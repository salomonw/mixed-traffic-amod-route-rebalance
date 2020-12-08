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


function add_demand_cnstr(model, G, g, xc, nodes, bush)
    links = keys(G)
    od_pairs = keys(g)
    if bush==0
        for j in nodes
            for od in od_pairs
                if od[1] == j 
                    @constraint(model, sum(xc[od][(i,j)] for i in get_incoming_edges(G, j)) + g[od] == sum(xc[od][(j,i)] for i in get_outgoing_edges(G, j)) )
                elseif od[2] == j
                    @constraint(model, sum(xc[od][(i,j)] for i in get_incoming_edges(G, j)) == g[od] + sum(xc[od][(j,i)] for i in get_outgoing_edges(G, j)) )
                else
                    @constraint(model, sum(xc[od][(i,j)] for i in get_incoming_edges(G, j)) == sum(xc[od][(j,i)] for i in get_outgoing_edges(G, j)) )
                end
            end
        end
    elseif bush==1
        O = unique([od[1] for od in od_pairs if g[od]>0])
        D = unique([od[2] for od in od_pairs if g[od]>0])
        od_nodes = unique([O;D])
        p = Dict()
        for j in od_nodes
            p[j] = sum(g[od] for od in od_pairs if od[2]==j) - sum(g[od] for od in od_pairs if od[1]==j)
        end
        # Global 
        for j in od_nodes
            @constraint(model, sum(xc[s][(i,j)] for s in O for i in get_incoming_edges(G,j) ) - sum( xc[s][(j,k)] for s in O for k in get_outgoing_edges(G,j) ) == p[j])
        end
        # Local 
        for s in O
            dsum = sum([g[od] for od in od_pairs if od[1]==s ])
            for j in nodes
                D = [od[2] for od in od_pairs if od[1]==s ]
                if j in D
                    d1 = g[(s,j)]
                    @constraint(model, sum([xc[s][(i,j)] for i in get_incoming_edges(G,j)]) - d1 == sum([xc[s][(j,k)] for k in get_outgoing_edges(G,j)]))
                elseif j == s
                    @constraint(model, sum([xc[s][(i,j)] for i in get_incoming_edges(G,j)]) == sum([xc[s][(j,k)] for k in get_outgoing_edges(G,j)])-dsum)
                else
                    @constraint(model, sum([xc[s][(i,j)] for i in get_incoming_edges(G,j)]) == sum([xc[s][(j,k)] for k in get_outgoing_edges(G,j)]))
                end
            end
        end
   end
end



#       for j in tnet.O:
#            p[j] = sum([tnet.g[(s,t)] for s,t in tnet.g.keys() if t==j]) - sum([tnet.g[(s,t)] for s,t in tnet.g.keys() if s==j])
        # Global
 #       for j in tnet.G_supergraph.nodes():
 #           m.addConstr(quicksum([x[(i,j)] for i,l in tnet.G_supergraph.in_edges(nbunch=j)]) - quicksum([x[(j,k)] for l,k in tnet.G_supergraph.out_edges(nbunch=j)]) == p[j] )
        # Local
  #      for s in tnet.O:
  #          dsum = sum([v for k,v in tnet.g.items() if k[0]==s])
  #          for j in tnet.G_supergraph.nodes():
  #              D = [d for o,d in tnet.g.keys() if o==s]
  #              if j in D:
  #                  d1 = tnet.g[(s,j)]
  #                  m.addConstr(quicksum(m.getVarByName('x^' + str(s) + '_' + str(i) + '_' + str(j)) for i, l in tnet.G_supergraph.in_edges(nbunch=j)) - d1 \
  #                      == quicksum(m.getVarByName('x^' + str(s) + '_' + str(j) + '_' + str(k)) for l, k in tnet.G_supergraph.out_edges(nbunch=j)))
  #              elif j == s:
  #                  m.addConstr(quicksum(m.getVarByName('x^' + str(s) + '_' + str(i) + '_' + str(j)) for i, l in tnet.G_supergraph.in_edges(nbunch=j)) \
  #                      == quicksum(m.getVarByName('x^' + str(s) + '_' + str(j) + '_' + str(k)) for l, k in tnet.G_supergraph.out_edges(nbunch=j)) - dsum)
  #              else:
  #                  m.addConstr(quicksum(m.getVarByName('x^' + str(s) + '_' + str(i) + '_' + str(j)) for i, l in tnet.G_supergraph.in_edges(nbunch=j)) \
  #                              == quicksum(m.getVarByName('x^' + str(s) + '_' + str(j) + '_' + str(k)) for l, k in tnet.G_supergraph.out_edges(nbunch=j)) )



function add_rebalancing_cnstr(m,G,xc,xr,nodes)
    od_pairs = keys(g)
    @constraint(m, [j in nodes], sum(xr[(i,j)] + sum(xc[od][(i,j)] for od in od_pairs) for i in get_incoming_edges(G, j)) == sum(xr[(j,k)] + sum(xc[od][(j,k)] for od in od_pairs) for k in get_outgoing_edges(G, j)))
end


function define_xc_vars(od_pairs, links, g, bush)
    xc = Dict()
    if bush==0
        for od in od_pairs
            xc[od] = Dict()
            for link in links
                xc[od][link] = @variable(m, lower_bound=0)
                set_name(xc[od][link], string("xc-(",od[1],",",od[2],")-(",link[1],",",link[2],")"))
            end
        end
    else
        O = unique([od[1] for od in od_pairs if g[od]>0])
        for o in O
            xc[o] = Dict()
            for link in links
                xc[o][link] = @variable(m, lower_bound=0)
                set_name(xc[o][link], string("xc-(",o,")-(",link[1],",",link[2],")"))
            end
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


function sol2dict(m, G, xc, od_pairs, links, bush)
    dict = Dict()
    for link in links
        flow = 0
        if bush == 0
            for od in od_pairs
                var = string("xc-(",od[1],",",od[2],")-(",link[1],",",link[2],")")
                flow +=  value(variable_by_name(m, var))
            end
        else
            for o in O
                var = string("xc-(",o,")-(",link[1],",",link[2],")")
                flow +=  value(variable_by_name(m, var))
            end
        end
        #reb_var = value(variable_by_name(m, string("xr-(",link[1],",",link[2],")")))
        dict[link] = Dict()
        dict[link]["flow"] = flow
       # dict[link]["flowNoRebalancing"] = flow - reb_var
   end
   return dict
end

function export_results(m, G, xc, od_pairs, links, fname, bush)
    dict = sol2dict(m, G, xc, od_pairs, links, bush)
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

bush = 1

m = Model(with_optimizer(Ipopt.Optimizer))

# Define variables
xc = define_xc_vars(od_pairs, links, g, bush)
#xr = define_xr_vars(links)


# Define sum od flow
O = unique([od[1] for od in od_pairs if g[od]>0])
x = Dict()
for link in links
    if bush==0
        x[link] = @NLexpression(m, sum(xc[od][link] for od in od_pairs))
    else
        x[link] = @NLexpression(m, sum(xc[o][link] for o in O))
    end
end

## Define Objective
@NLobjective(m, Min, sum( x[link]* G[link]["t_0"]* sum( fcoeffs[n]*((x[link] + G_exogenous[link]["flow"])/G[link]["capacity"])^(n-1) for n=1:N ) for link in links))

## Define Constraints
add_demand_cnstr(m, G, g, xc, nodes, bush)

optimize!(m)

export_results(m, G, xc, od_pairs, links, out_file, bush)

# write runnning time to file
solvetime = MOI.get(m, MOI.SolveTime())
ofile= open("tmp/solvetime.txt", "w")
println(ofile, solvetime)
close(ofile)
