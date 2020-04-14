using JSON
using JuMP, Ipopt
using Base

function read_json_file(fname)
    #TODO: write description
    dict2 = Dict()
    open(fname, "r") do f
        dicttxt = read(f, String)
        dict2 = JSON.parse(dicttxt)
    end
    return dict2
end

function read_nx_file(fname, flow)
    #TODO: write description
    G_dict = read_json_file(fname)
    G = Dict()
    nodes = []
    for edge in G_dict["links"]
        s = string(edge["source"])
        t = string(edge["target"])
        G[(s,t)] = Dict()
        push!(nodes, s)
        push!(nodes, t)
        nodes = unique(nodes)

        G[(s,t)]["source"] = s
        G[(s,t)]["target"] = t

        if flow == 0
            for att in ["capacity", "length", "t_0"]
                G[(s,t)][att] = edge[att]

            end
        else
            for att in ["capacity", "length", "t_0", "flow"]
                G[(s,t)][att] = edge[att]
            end            
        end
    end
    return G, sort(nodes)

end


function read_fcoeffs(fname)
    #TODO: write description
    return read_json_file(fname)
end



function read_g(fname)
    #TODO: write description
    g_dict = Dict()
    dict2 = read_json_file(fname)
    for od in keys(dict2)
        a = split(split(split(od, "(")[2],")")[1],",")
        g_dict[(a[1], a[2][2:end])] = dict2[od]
    end
    return g_dict
end


function dict2json(dict, fname)
    json_string = JSON.json(dict)
    open(fname,"w") do f 
        write(f, json_string) 
    end
end


function get_incoming_edges(G, i)
    a = []
    a = unique([if G[k]["target"]==i G[k]["source"] else 0 end for k in keys(G)])
    filter!(x->x≠0,a)
    return a
end



function get_outgoing_edges(G, i)
    a = []
    a = unique([if G[k]["source"]==i G[k]["target"] else 0 end for k in keys(G)])
    filter!(x->x≠0,a)
    return a
end


function add_demand_cnstr(model,G,g,xc,nodes)
    od_pairs = keys(g)
    links = keys(G)
    for j in nodes
        for od in od_pairs
            if od[1] == j 
                @constraint(m, sum(xc[od][(i,j)] for i in get_incoming_edges(G, j)) + g[od] == sum(xc[od][(j,i)] for i in get_outgoing_edges(G, j)) )
            elseif od[2] == j
                @constraint(m, sum(xc[od][(i,j)] for i in get_incoming_edges(G, j)) == g[od] + sum(xc[od][(j,i)] for i in get_outgoing_edges(G, j)) )
            else
                @constraint(m, sum(xc[od][(i,j)] for i in get_incoming_edges(G, j)) == sum(xc[od][(j,i)] for i in get_outgoing_edges(G, j)) )
            end
        end
    end
end


function add_rebalancing_cnstr(m,G,xc,xr,nodes)
    od_pairs = keys(g)
    @constraint(m, [j in nodes], sum(xr[(i,j)] + sum(xc[od][(i,j)] for od in od_pairs) for i in get_incoming_edges(G, j)) == sum(xr[(j,k)] + sum(xc[od][(j,k)] for od in od_pairs) for k in get_outgoing_edges(G, j)))
end


function define_xc_vars(od_pairs, links)
    xc = Dict()
    for od in od_pairs
        xc[od] = Dict()
        for link in links
            xc[od][link] = @variable(m, lower_bound=0)
            set_name(xc[od][link], string("xc-(",od[1],",",od[2],")-(",link[1],",",link[2],")"))
        end
    end
    return xc
end

function define_xr_vars(links)
    xr = Dict()
    for link in links
        xr[link] = @variable(m, lower_bound=0)
        set_name(xr[link], string("xr-(",link[1],",",link[2],")"))
    end
    return xr
end


function sol2dict(m, G, xc, od_pairs, links)
    dict = Dict()
    for link in links
        flow = 0
        for od in od_pairs
            var = string("xc-(",od[1],",",od[2],")-(",link[1],",",link[2],")")
            flow +=  value(variable_by_name(m, var))
        end
        #reb_var = value(variable_by_name(m, string("xr-(",link[1],",",link[2],")")))
        dict[link] = Dict()
        dict[link]["flow"] = flow
       # dict[link]["flowNoRebalancing"] = flow - reb_var
    end
    return dict
end

function export_results(m, G, xc, od_pairs, links, fname)
    dict = sol2dict(m, G, xc, od_pairs, links)
    dict2json(dict, fname)
end

function print_sol(m, xc, od_pairs, links, flag=1)
    if flag == 1
        for od in od_pairs
            xc[od] = Dict()
            for link in links
                var = string("xc-(",od[1],",",od[2],")-(",link[1],",",link[2],")")
                println(var, "= ", value(variable_by_name(m, var)))
            end
        end
        
        for link in links
            var = 
            println(value(variable_by_name(m, string("xr-(",link[1],",",link[2],")"))))
        end
    else
        return 0
    end
end





#G, nodes = read_nx_file("tmp/G.json", 0)
G, nodes = read_nx_file("tmp/G_supergraph.json", 0)
G_exogenous, nodes_ex = read_nx_file("tmp/exogenous_G.json", 1)
#G_exogenous = "False"
fcoeffs = read_fcoeffs("tmp/fcoeffs.json")
g = read_g("tmp/g.json")
out_file = "tmp/out.json"


od_pairs = keys(g)
links = keys(G)
N = length(fcoeffs)

m = Model(with_optimizer(Ipopt.Optimizer))

# Define variables
xc = define_xc_vars(od_pairs, links)
#xr = define_xr_vars(links)

# Define sum od flow
x = Dict()
for link in links
    x[link] = @NLexpression(m, sum(xc[od][link] for od in od_pairs))
end

## Define Objective
@NLobjective(m, Min, sum( x[link]* G[link]["t_0"]* sum( fcoeffs[n]*((x[link] + G_exogenous[link]["flow"])/G[link]["capacity"])^(n-1) for n=1:N ) for link in links))

## Define Constraints
add_demand_cnstr(m, G, g, xc, nodes)
#add_rebalancing_cnstr(m, G, xc, xr, nodes)

optimize!(m)

export_results(m, G, xc, od_pairs, links, out_file)

# write runnning time to file
solvetime = MOI.get(m, MOI.SolveTime())
ofile= open("tmp/solvetime.txt", "w")
println(ofile, solvetime)
close(ofile)
