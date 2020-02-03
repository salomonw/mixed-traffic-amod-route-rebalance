
import networkx as nx
import scipy.optimize

def BPR(x, m, fcoeffs, exo):
	return sum([fcoeffs[i]*((x+exo)/m)**i for i in range(len(fcoeffs))])


def BPR_social(x, m, fcoeffs, exo):
	return sum([fcoeffs[i]*((x+exo)/m)**i for i in range(len(fcoeffs))]) + sum([fcoeffs[i+1]*(i+1) *((x+exo)/m)**(i) for i in range(len(fcoeffs)-1)])

def set_zero_flow_x0(G):
	for (i,j) in G.edges():
		G[i][j]['flow'] = 0
		G[i][j]['x'] = 0
		G[i][j]['t_k'] = G[i][j]['t_0']

def set_zero_x0(G):
	for (i,j) in G.edges():
		G[i][j]['x'] = 0
		G[i][j]['t_k'] = G[i][j]['t_0']

def all_or_nothing(G, g, fcoeffs, a, G_exo=False):
	for (i,j) in G.edges():
		G[i][j]['flow'] = a * G[i][j]['x'] + (1-a) * G[i][j]['flow']
		# update travel times
		if G_exo:
			exo = G_exo[i][j]['flow']
		else:
			exo=0
		G[i][j]['t_k'] = G[i][j]['t_0'] * BPR( G[i][j]['flow'] ,G[i][j]['capacity'], fcoeffs, exo)
		G[i][j]['x'] = 0

def all_or_nothing_social(G, g, fcoeffs, a, G_exo=False):
	for (i,j) in G.edges():
		G[i][j]['flow'] = a * G[i][j]['x'] + (1-a) * G[i][j]['flow']
		# update travel times
		if G_exo:
			exo = G_exo[i][j]['flow']
		else:
			exo=0
		G[i][j]['t_k'] = G[i][j]['t_0'] * BPR_social( G[i][j]['flow'] ,G[i][j]['capacity'], fcoeffs, exo)
		G[i][j]['x'] = 0

def update_stepsize_MSA(iteration):
	a = 1/((1+iteration))
	return a

def get_total_system_travel_time(G, fcoeffs):
	return sum([G[i][j]['flow'] * G[i][j]['t_k'] for (i,j) in G.edges() ])

def get_shortest_path_travel_time(G, g, fcoeffs):
	sptt = 0
	for (s,t) in g.keys():
		path = nx.shortest_path(G, source=s, target=t, weight='t_k')
		for i, j in zip(path, path[1:]):
			sptt += G[i][j]['t_k']*g[(s,t)]
			G[i][j]['x'] += g[(s,t)]
	return sptt

def update_stepsize_FW(G, fcoeffs, j):
	#TODO: add the G_exo part
	sol = scipy.optimize.root(derivative_FW, x0=1/(j+1), args=(G, fcoeffs), jac=False)
	a = max(0, min(1, sol.x[0]))
	return a

def derivative_FW(a, G, fcoeffs, G_exo=False):
	sum_derivative = 0
	for i,j in G.edges():
		if G_exo:
			exo=G_exo[i][j]['flow']
		else:
			exo = 0
		sum_derivative += (G[i][j]['flow']-G[i][j]['x']) * G[i][j]['t_0'] *\
		BPR( G[i][j]['flow'] + a*(G[i][j]['flow']-G[i][j]['x']) , G[i][j]['capacity'], fcoeffs, exo )
	return sum_derivative


def assignment(G, g, fcoeffs, G_exo=False, flow=True, method='FW', accuracy=0.0001, max_iter=2000):
	if flow==True:
		set_zero_x0(G)
	else:
		set_zero_flow_x0(G)
	RG = 100000
	j = 0
	while RG > accuracy and j<max_iter:	
		sptt = get_shortest_path_travel_time(G,g,fcoeffs)
		if method=='MSA':
			a = update_stepsize_MSA(j)
		if method=='FW':
			a = update_stepsize_FW(G, fcoeffs, j)
		all_or_nothing(G, g, fcoeffs, a, G_exo)
		tstt = get_total_system_travel_time(G, fcoeffs)
		RG = abs(tstt/sptt - 1)
		j += 1
		#print(RG)
	return G


def assignment_social(G, g, fcoeffs, G_exo=False, flow=True, method='FW', accuracy=0.0001, max_iter=1000):
	G1 = G.copy
	if flow==True:
		set_zero_x0(G1)
	else:
		set_zero_flow_x0(G1)
	RG = 100000
	j = 0
	while RG > accuracy and j<max_iter:
		sptt = get_shortest_path_travel_time(G1,g,fcoeffs)
		if method=='MSA':
			a = update_stepsize_MSA(j)
		if method=='FW':
			a = update_stepsize_FW(G1, fcoeffs, j)
		all_or_nothing_social(G1, g, fcoeffs, a, G_exo)
		tstt = get_total_system_travel_time(G1, fcoeffs)
		RG = abs(tstt/sptt - 1)
		j += 1
		#print(RG)
	return G1