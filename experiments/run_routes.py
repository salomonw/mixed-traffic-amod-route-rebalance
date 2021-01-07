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
    else:
        netFile, gFile, fcoeffs, tstamp, dir_out = tnet.get_network_parameters(net_name=net_name,
                                                                               experiment_name=net_name + '_n_variation')
        tNet = tnet.tNet(netFile=netFile, gFile=gFile, fcoeffs=fcoeffs)
    return tNet, fcoeffs


def od_travel_times(tNet, s_flows, od):
	routes = routeFinder.RouteFinder(tNet.G, tNet.g, s_flows, eps=0.0, od=od)
	rs = {}
	for i, route in routes.items():
		rs[i]={'p': route['p'], 'tt':routeFinder.RouteTravelTime(tNet.G,route['r']), 'r':route['r']}
	return rs



def plot_network(G, ax, width=1, cmap=plt.cm.Blues, edge_width=False, 
	edgecolors=False, nodecolors=False, nodesize=False, arrowsize=False):
    pos = nx.get_node_attributes(G, 'pos')
    nx.draw(G, pos, width=edge_width,  ax=ax, edge_color=edgecolors, node_size=nodesize, node_color=nodecolors,
    	connectionstyle='arc3, rad=0.04',arrowsize=2, arrowstyle='fancy')



def plot_routes(G, od, result, ax):
	#fig, ax = plt.subplots()
	colors = ['r', 'b', 'y', 'orange']
	tNet.read_node_coordinates('data/pos/'+net+'.txt')
	#tnet.plot_network(tNet.G, width=0.3)
	node_color = ['b' if n in [od[0]] else 'gray' for n in G.nodes()]
	node_color = ['g' if n in [od[1]] else 'gray' for n in G.nodes()]
	node_size = [25 if n in [od[0], od[1]] else 5 for n in G.nodes()]
	l =0
	for i, dic in result.items():
		r = dic['r']
		r_links = [(r[n],r[n+1]) for n in range(len(r)-1)]
		#for link in r_links:
		edge_colors = [colors[l] if e in r_links else 'gray' for e in G.edges()]
		edge_width = [2 if e in r_links else 0.2 for e in G.edges()]
		arrow_size = [3 if e in r_links else 0.2 for e in G.edges()]
		plot_network(tNet.G, ax, edge_width=edge_width, 
			edgecolors=edge_colors, nodecolors=node_color, 
			nodesize=node_size, arrowsize=arrow_size)
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



#netname = 'NYC'
net = 'EMA_mid'
g_mul = 2

tNet, fcoeffs = read_net(net)
#tNet_UC = copy.deepcopy(tNet)

tNet.build_supergraph(identical_G=True)
g_per = tnet.perturbDemandConstant(tNet.g, g_mul)
tNet.set_g(g_per)


#tNet_UC.build_supergraph(identical_G=True)
#g_per = tnet.perturbDemandConstant(tNet_UC.g, g_mul)
#tNet_UC.set_g(g_per)

tNet, runtime, s_flows = cars.solve_bush_CARSn(tNet, fcoeffs=fcoeffs, n=8,
                                                exogenous_G=False, rebalancing=False,
                                                linear=False, bush=True)

#tNet_UC, runtime, s_flows_UC = cars.solve_bush_CARSn(tNet_UC, fcoeffs=fcoeffs, n=8,
#                                                exogenous_G=False, rebalancing=False,
#                                                linear=False, bush=True, userCentric=True)


cars.supergraph2G(tNet)
#cars.supergraph2G(tNet_UC)

objSO = sum([tNet.G[i][j]['flow']*tNet.G[i][j]['t_k'] for i,j in tNet.G.edges()])
#objUC = sum([tNet_UC.G[i][j]['flow']*tNet_UC.G[i][j]['t_k'] for i,j in tNet_UC.G.edges()])

print(objSO)
#print(objUC)
#print(objUC/objSO)
#rebRoutes = routeFinder.rebRouteFinder(tNet.G, eps=0)
#userRoutes = routeFinder.userRouteFinder(tNet.G, tNet.g, s_flows, eps=0)

#select OD pair
random.seed(8)
ods = dict(sorted(tNet.g.items(), key=lambda item: item[1]))
ods = list(ods.keys())[-30:]

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

for od in ods[0:10]:
	result = od_travel_times(tNet, s_flows, od)
	fig, ax = plt.subplots()
	plot_routes(tNet.G_supergraph, od, result, ax)
	plt.savefig('plot_'+str(od)+'.pdf')	
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


