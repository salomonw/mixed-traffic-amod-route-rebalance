


def run_TA(G, demand, fcoeffs):
    net = nx_2_net(G)
    



def nx2vec(G):
    link_id = {}
    idx = 0
    for (s,t) in G.edges():        
        link_id[idx] = (s,t)
        net[0][idx] = s
        net[1][idx] = t
        net[2][idx] = G[s][t]['capacity']
        net[3][idx] = G[s][t]['length']
        net[4][idx] = G[s][t]['t_0']
        net[5][idx] = G[s][t]['B']
        net[6][idx] = G[s][t]['power']
        net[7][idx] = G[s][t]['speedLimit']
        net[8][idx] = G[s][t]['toll']
        net[9][idx] = G[s][t]['type']
    return net, link_id



#y = CC.shortest_path(net,od, [0]*76)

#print "this", y.label_correcting(1)[3]
#print "links visited", y.visited_links(1, 21, y.label_correcting(1)[3])
start_time = time.time()

z = CC.stochastic_loading(net, od, [0]*76)
#for a, b in od.items():
#print z.dial(2)

a = CC.assignment(net, od, [0]*76)
#print a.msa_ue()
print(len(net[0]))
print("--- %s seconds ---" % (time.time() - start_time))


