

using JSON
using JuMP, Ipopt

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
    nodes = Int64[]
    for edge in G_dict["links"]
        G[(edge["source"],edge["target"])] = Dict()
        push!(nodes, edge["source"])
        push!(nodes, edge["target"])
        nodes = unique(nodes)

        if flow == 0
            for att in ["capacity", "length", "t_0"]
                G[(edge["source"],edge["target"])][att] = edge[att]

            end
        else
            for att in ["capacity", "length", "t_0", "flow"]
                G[(edge["source"],edge["target"])][att] = edge[att]
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
        g_dict[(parse(Int64, a[1]), parse(Int64, a[2]))] = dict2[od]
    end
    return g_dict
end


function get_incoming_edges(G, i)

end



function get_outgoing_edges(G, i)

end

G, nodes = read_nx_file("tmp/G.json", 0)
G_exogenous = read_nx_file("tmp/exogenous_G.json", 1)
fcoeffs = read_fcoeffs("tmp/fcoeffs.json")
g = read_g("tmp/g.json")

print(nodes)
od_pairs = keys(g)
links = keys(G)
N = length(fcoeffs)

numLinks = length(G)

m = Model(with_optimizer(Ipopt.Optimizer))

# Define variables
xc = Dict()
for od in od_pairs
    xc[od] = Dict()
    for link in links
        xc[od][link] = @variable(m)#base_name=string("xc-",parse(od,String), "-", parse(link,String)))
    end
end


# Define Objective
@NLobjective(m, Min, sum(  sum(xc[od][link] for od in od_pairs) * G[link]["t_0"] * sum( fcoeffs[n]* sum(xc[od][link] for od in od_pairs)^(n-1) for n=1:N ) for link in links))

# Define Constraints
for j in nodes
    for od in od_pairs
        if od[1] == j
            @constraint(x[od][j])
        elseif od[2] == j

        else
    end
end


for j in tnet.G_supergraph.nodes():
    for w, d in tnet.g.items():
        if j == w[0]:
            m.addConstr(quicksum(m.getVarByName('x^'+str(w)+'_'+str(i)+'_'+str(j)) for i,l in tnet.G_supergraph.in_edges(nbunch=j)) + d == quicksum(m.getVarByName('x^'+str(w)+'_'+str(j)+'_'+str(k)) for l,k in tnet.G_supergraph.out_edges(nbunch=j)))
        elif j == w[1]:
            m.addConstr(quicksum(m.getVarByName('x^' + str(w) + '_' + str(i) + '_' + str(j)) for i,l in tnet.G_supergraph.in_edges(nbunch=j)) == quicksum(m.getVarByName('x^' + str(w) + '_' + str(j) + '_' + str(k)) for l,k in tnet.G_supergraph.out_edges(nbunch=j)) + d)
        else:
            m.addConstr(quicksum(m.getVarByName('x^' + str(w) + '_' + str(i) + '_' + str(j)) for i,l in tnet.G_supergraph.in_edges(nbunch=j)) == quicksum(m.getVarByName('x^' + str(w) + '_' + str(j) + '_' + str(k)) for l,k in tnet.G_supergraph.out_edges(nbunch=j)))

m.update()


#print(obj)
#@NLobjective()
#@NLobjective(m, Min, obj)

# Define constraints
#@constraint(m::Model, expr)