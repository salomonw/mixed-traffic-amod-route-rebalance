import src.tnet as tnet
import src.CARS as cars
import src.routeFinder as routeFinder
import experiments.build_NYC_subway_net as nyc
import random
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import math 
import pandas as pd
import copy 

def read_net(net_name):
    if net_name == 'NYC':
        tNet, tstamp, fcoeffs = nyc.build_NYC_net('data/net/NYC/', only_road=True)
        tNet.build_supergraph(identical_G=True)
        tNet.read_node_coordinates('data/pos/NYC.txt')
    else:
        netFile, gFile, fcoeffs, tstamp, dir_out = tnet.get_network_parameters(net_name=net_name,
                                                                               experiment_name=net_name + '_n_variation')
        tNet = tnet.tNet(netFile=netFile, gFile=gFile, fcoeffs=fcoeffs)
    return tNet, fcoeffs


def od_travel_times(tNet, s_flows, od):
	routes = routeFinder.RouteFinder(tNet.G_supergraph, tNet.g, s_flows, eps=10.0, od=od)
	rs = {}
	for i, route in routes.items():
		rs[i]={'p': route['p'], 'tt':routeFinder.RouteTravelTime(tNet.G_supergraph,route['r']), 'r':route['r']}
	return rs



def plot_network(G, ax, width=1, cmap=plt.cm.Blues, edge_width=False, 
	edgecolors=False, nodecolors=False, nodesize=False, arrowsize=False,edge_alpha=1,linkstyle='-' ):
    pos = nx.get_node_attributes(G, 'pos')
    if linkstyle == '-':
        nx.draw(G, pos, 
                width=edge_width,  
                ax=ax, 
                edge_color=edgecolors, 
                node_size=nodesize, 
                node_color=nodecolors,
                #connectionstyle='arc3, rad=0.04',
                arrowsize=0.5, 
                arrowstyle='fancy',
                alpha=edge_alpha)
    else:
        nx.draw(G, pos, width=edge_width,  
                ax=ax, edge_color=edgecolors, 
                node_size=nodesize, node_color=nodecolors,
                #connectionstyle='arc3, rad=0.04',
                arrowsize=0.5, arrowstyle='fancy',
                alpha=edge_alpha, style=linkstyle)


def plot_routes(tNet, od, result, ax):
	#fig, ax = plt.subplots()
	cmap2 =  ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])
             for i in range(30)]
	cmap = ['b','m','y','c','g', 'r']  
	cmap.extend(cmap2)
	tNet.read_node_coordinates('data/pos/'+net+'.txt')
	#tnet.plot_network(tNet.G, width=0.3)
	node_color = ['red' if n in [od[0]] else ('green' if n in [od[1]] else 'gray') for n in tNet.G.nodes()]
	#node_color = ['green' if n in [od[1]] else 'gray' for n in tNet.G_supergraph.nodes()]
	node_size = [25 if n in [od[0], od[1]] else 0.2 for n in tNet.G.nodes()]
	l =0
	plot_network(tNet.G, ax, edge_width=0.35,
                     edgecolors='gray', nodecolors=node_color,
                     nodesize=node_size, arrowsize=0.15,edge_alpha=0.7)
    
	for i, dic in result.items():
		print(i)
		r = dic['r']
		r_links = [(r[n],r[n+1]) for n in range(len(r)-1)]
		G_ = tNet.G_supergraph.edge_subgraph(r_links)
        
		edges = ((edge[0],edge[1]) for edge in G_.edges(data=True) if edge[2]['type']=="'")	
		G = G_.edge_subgraph(edges)
		edge_colors = [cmap[l] if e in r_links else 'gray' for e in G.edges()]
		edge_width = [2 if e in r_links else 0 for e in G.edges()]
		arrow_size = [1 if e in r_links else 0 for e in G.edges()]
		plot_network(G, ax, edge_width=edge_width, 
			edgecolors=edge_colors, nodecolors='gray', 
			nodesize=0.0, arrowsize=arrow_size,edge_alpha=0.7, 
            linkstyle='-.')

		edges = ((edge[0],edge[1]) for edge in G_.edges(data=True) if edge[2]['type']=="s")
		G = G_.edge_subgraph(edges)
		edge_colors = [cmap[l] if e in r_links else 'gray' for e in G.edges()]
		edge_width = [1 if e in r_links else 0 for e in G.edges()]
		arrow_size = [1 if e in r_links else 0 for e in G.edges()]
		plot_network(G, ax, edge_width=edge_width, 
			edgecolors=edge_colors, nodecolors='gray', 
			nodesize=0.0, arrowsize=arrow_size,edge_alpha=0.7, 
            linkstyle='--')        
        
		edges = ((edge[0],edge[1]) for edge in G_.edges(data=True) if edge[2]['type']=="b")
		G = G_.edge_subgraph(edges)
		edge_colors = [cmap[l] if e in r_links else 'gray' for e in G.edges()]
		edge_width = [0.5 if e in r_links else 0 for e in G.edges()]
		arrow_size = [1 if e in r_links else 0 for e in G.edges()]
		plot_network(G, ax, edge_width=edge_width, 
			edgecolors=edge_colors, nodecolors='gray', 
			nodesize=0.1, arrowsize=arrow_size,edge_alpha=0.7, 
			linkstyle=(0, (1, 10)))
        
		edges = ((edge[0],edge[1]) for edge in G_.edges(data=True) if edge[2]['type']=="p")
		G = G_.edge_subgraph(edges)
		edge_colors = [cmap[l] if e in r_links else 'gray' for e in G.edges()]
		edge_width = [2 if e in r_links else 0 for e in G.edges()]
		arrow_size = [1 if e in r_links else 0 for e in G.edges()]
		plot_network(G, ax, edge_width=edge_width, 
			edgecolors=edge_colors, nodecolors='gray', 
			nodesize=0.0, arrowsize=arrow_size,edge_alpha=0.7, 
            linkstyle='-.')
        
		edges = ((edge[0],edge[1]) for edge in G_.edges(data=True) if edge[2]['type']==0)
		G = G_.edge_subgraph(edges)
		edge_colors = [cmap[l] if e in r_links else 'gray' for e in G.edges()]
		edge_width = [2 if e in r_links else 0 for e in G.edges()]
		arrow_size = [1 if e in r_links else 0 for e in G.edges()]
		plot_network(G, ax, edge_width=edge_width, 
			edgecolors=edge_colors, nodecolors='gray', 
			nodesize=0.0, arrowsize=arrow_size,edge_alpha=0.7, 
            linkstyle='-')
		l+=1
	return ax

def weighted_avg_and_std(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, weights=weights)
    # Fast and numerically precise:
    variance = np.average((values-average)**2, weights=weights)
    return math.sqrt(variance)


#net = 'EMA_mid'
net = 'NYC'
g_mul = 1.5

# add road network
tNet, fcoeffs = read_net(net)
# add biking layer
tNet.build_layer(one_way=True, avg_speed=6,capacity=99999, symb="b", identical_G=False)
# add pedestrian network
tNet.build_layer(one_way=False, avg_speed=3.1, capacity=99999, symb="'", identical_G=False)
# add subway for nyc
layer = tnet.readNetFile(netFile='data/net/NYC/NYC_M_Subway_net.txt')
tNet.add_layer(layer=layer, layer_symb='s')

# set demand
g_per = tnet.perturbDemandConstant(tNet.g, g_mul)
tNet.set_g(g_per)

# solve CARS model
tNet, runtime, s_flows = cars.solve_bush_CARSn(tNet, fcoeffs=fcoeffs, n=8,
                                                exogenous_G=False, rebalancing=False,
                                                linear=False, bush=True)


#select OD pair
random.seed(20)
ods = list(tNet.g.keys())
random.shuffle(ods)
ods = list(ods)[-50:]

#ods = [(269,546), (1034,899), (1034,958)]

table = {}
table['od'] = []
table['d'] = []
table['n'] = []
table['n_UC'] = []
table['spTT'] = []
table['avgTT'] = []
table['WavgTT'] = []
table['TT_UC'] = []
table['stdTT'] = []
table['WstdTT'] = []

for od in ods:
	result = od_travel_times(tNet, s_flows, od)
	#print(result)
	fig, ax = plt.subplots()
	plot_routes(tNet, od, result, ax)
	#plt.show()
	plt.savefig('plot_'+net+'_'+str(od)+'.pdf')	
	'''
	resultUC = od_travel_times(tNet_UC, s_flows_UC, od)
	print([v['tt'] for k,v in resultUC.items()])
	tts = [(v['p'], v['tt']) for k,v in result.items()]
	ps = [v['p'] for k,v in result.items()]
	sp = nx.shortest_path(tNet.G, od[0],od[1], weight='t_k')
	sp_UC = nx.shortest_path(tNet_UC.G, od[0],od[1], weight='t_k')
	table['od'].append(od)
	table['d'].append(tNet.g[od])
	table['n'].append(len(result))
	table['n_UC'].append(len(resultUC))
	table['avgTT'].append(np.average(tts))
	table['spTT'].append(routeFinder.RouteTravelTime(tNet.G,sp))
	table['WavgTT'].append(np.average(tts, weights=ps))
	table['TT_UC'].append(routeFinder.RouteTravelTime(tNet_UC.G,sp_UC))
	table['stdTT'].append(np.std(tts))
	table['WstdTT'].append(weighted_avg_and_std(tts, ps))
	df = pd.DataFrame.from_dict(table, orient='index').transpose()
	print(df)
	'''


