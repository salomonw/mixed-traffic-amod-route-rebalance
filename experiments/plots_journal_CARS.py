import src.tnet as tnet
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import math    

plt.style.use(['science','ieee', 'high-vis'])


def txt2list(fname):
	return [line for line in open(fname)]

def read_result(fname):
	df = pd.read_csv(fname)
	results = df.T.values.tolist()
	return results

def read_parameters(fname):
	dic = {}
	for line in open(fname, 'r').readlines():
		p,v = line.split()
		dic[p] = v
	return dic


def plot_topology(netname):
    netFile, gFile, fcoeffs, tstamp, dir_out = tnet.get_network_parameters(net_name=netname,

                                                                           experiment_name=netname + 'topo_plot')
    tNet = tnet.tNet(netFile=netFile, gFile=gFile, fcoeffs=fcoeffs)
    tNet.read_node_coordinates('data/pos/'+netname+'.txt')
    fig, ax = tnet.plot_network(tNet.G, width=0.3)
    return fig, ax


def plot_convergance(fname_sys, fname_usr):
	return 1


def plot_costPenRate(fname, ax, parameters, k):
	j, cavsCost, noncavsCost, totCost, cavsFlow, nonCavsFlow, pedestrianFlow, rebalancingFlow, bikeFlow, subwayFlow = read_result(fname)
	if k == 'A':
		for i in range(len(cavsCost)):
			cavsCost[i] = max(noncavsCost[i], cavsCost[i])
			totCost[i] = max(noncavsCost[i], totCost[i])
	j = [round(.1 * i, 1) for i in range(11)]
	lstyle = ['-', '--', ':']
	i = 0
	alg = 'CARS'+parameters['n:']
	ax.plot(j, noncavsCost, label='Private', linestyle=lstyle[i], linewidth=2, marker='x')
	ax.plot(j, cavsCost, label='AMoDs', linestyle=lstyle[i], linewidth=2, marker="^")
	ax.plot(j, totCost, label='Total', linestyle=lstyle[i], linewidth=2, marker='o')
	ax.legend()
	ax.set_xlabel('Penetration Rate')
	ax.set_ylabel('Avg. Travel Time (min)')
	ax.set_xlim((0, 1))
	ax.legend(framealpha=0.8, fontsize='small', frameon=True, facecolor='w', fancybox='False')
	#ax.legend.get_frame().set_linewidth(0.2)
	return ax


def plot_flowPenRate(fname, ax, parameters):
	n, cavsCost, noncavsCost, totCost, cavsFlow, nonCavsFlow, pedestrianFlow, rebalancingFlow, bikeFlow, subwayFlow = read_result(fname)
	width = 0.9
	x_name = [round(.1 * i, 1) for i in range(11)]
	x = list(range(len(x_name)))
	p1 = ax.bar(x, nonCavsFlow, width, label='Private')
	p2 = ax.bar(x, cavsFlow, width,
	             bottom=nonCavsFlow, label='AMoD')
	p3 = ax.bar(x, rebalancingFlow, width,
	             bottom=[cavsFlow[i] + nonCavsFlow[i] for i in range(len(cavsFlow))], label='Rebalancing')
	if sum(subwayFlow)>10:
		p6 = ax.bar(x, subwayFlow, width,
	             bottom=[cavsFlow[i] + nonCavsFlow[i] + rebalancingFlow[i] + pedestrianFlow[i] + bikeFlow[i] for i in
	                     range(len(cavsFlow))], label='Subway')
	if sum(pedestrianFlow)>10:
		p4 = ax.bar(x, pedestrianFlow, width,
	             bottom=[cavsFlow[i] + nonCavsFlow[i] + rebalancingFlow[i] for i in range(len(cavsFlow))], label='Pedestrian')
	if sum(bikeFlow)>10:
		p5 = ax.bar(x, bikeFlow, width,
	             bottom=[cavsFlow[i] + nonCavsFlow[i] + rebalancingFlow[i] + pedestrianFlow[i] for i in
	                     range(len(cavsFlow))], label='Biking')


	ax.set_ylabel('Miles per mode of transport')
	ax.set_xlabel('Penetration rate')
	ax.set_xticks(x)
	ax.set_xticklabels(x_name)
	ax.legend(framealpha=0.8, fontsize='small', frameon=True, loc=3, facecolor='w', fancybox='False')
	#ax.legend.get_frame().set_linewidth(0.2)
	return ax

'''
dire = '2021-01-08_11:51:44_penRate_NYC_1.5ASB_Reb_True'
fname = 'results/' + dire + '/results.csv' 
parameters = read_parameters('results/' + dire + '/parameters.txt' )
#print(read_result(fname))

fig, ax = plt.subplots(1 ,figsize=(2.5,2))
plot_costPenRate(fname, ax, parameters)
plt.savefig('a.pdf')

fig, ax = plt.subplots(1 ,figsize=(3.6,2))
plot_flowPenRate(fname, ax, parameters)
plt.savefig('b.pdf')
'''

# comparison

def plot_comparison(fnames, out):
	fig, ax = plt.subplots(ncols=2, 
							nrows=len(fnames), 
						#	width_ratios=[1,2], 
							gridspec_kw={'width_ratios':[1,2]},
							figsize=(3.6*2, 2*len(fnames)), 
							#sharex=True, 
							sharey=False)
	j = 0
	for f in fnames:
		fname = 'results/' + f + '/results.csv'
		parameters = read_parameters('results/' + f + '/parameters.txt' )
		if out =='1c':
			plot_costPenRate(fname, ax[j,0], parameters, 'A')
		else:
			plot_costPenRate(fname, ax[j,0], parameters, 'B')
		plot_flowPenRate(fname, ax[j,1], parameters)
		j  +=1
	#plt.legend(frameon=True, fancybox=False)
	plt.tight_layout()
	plt.savefig(out+'.pdf')
	#plt.show()

one = '2021-01-08_11/50/19_penRate_NYC_1.0A_Reb_True'.replace('/', ':')
two = '2021-01-08_11/50/08_penRate_NYC_1.5A_Reb_True'.replace('/', ':')
three = '2021-01-08_11/51/44_penRate_NYC_2.0A_Reb_True'.replace('/', ':')
four = '2021-01-08_11/51/44_penRate_NYC_4.0A_Reb_True'.replace('/', ':')
fnames = [one, two, three, four]

plot_comparison(fnames,'1c')


one = '2021-01-08_11/50/19_penRate_NYC_1.0AS_Reb_True'.replace('/', ':')
two = '2021-01-08_11/50/08_penRate_NYC_1.5AS_Reb_True'.replace('/', ':')
three = '2021-01-08_11/51/44_penRate_NYC_2.0AS_Reb_True'.replace('/', ':')
four = '2021-01-08_11/51/43_penRate_NYC_4.0AS_Reb_True'.replace('/', ':')
fnames = [one, two, three, four]

plot_comparison(fnames,'1_5c')




one = '2021-01-08_11/50/08_penRate_NYC_1.0ASP_Reb_True'.replace('/', ':')
two = '2021-01-08_11/51/48_penRate_NYC_1.5ASP_Reb_True'.replace('/', ':')	
three = '2021-01-08_11/51/44_penRate_NYC_2.0ASP_Reb_True'.replace('/', ':')
four = '2021-01-08_11/52/40_penRate_NYC_4.0ASP_Reb_True'.replace('/', ':')
fnames = [one, two, three, four]

plot_comparison(fnames,'2c')



one = '2021-01-08_11/50/08_penRate_NYC_1.0ASPB_Reb_True'.replace('/', ':')
two = '2021-01-08_11/51/44_penRate_NYC_1.5ASPB_Reb_True'.replace('/', ':')
three = '2021-01-12_00:58:41_penRate_NYC_2.0ASPB_Reb_True'.replace('/', ':')
four = '2021-01-14_02:00:28_penRate_NYC_4.0ASPB_Reb_True'.replace('/', ':')
fnames = [one, two, three, four]

plot_comparison(fnames,'4c')

one = '2021-01-08_11/51/44_penRate_NYC_2.0A_Reb_True'.replace('/', ':')
two = '2021-01-08_11/51/44_penRate_NYC_2.0AS_Reb_True'.replace('/', ':')
three = '2021-01-08_11/51/44_penRate_NYC_2.0ASP_Reb_True'.replace('/', ':')
four = '2021-01-12_00:58:41_penRate_NYC_2.0ASPB_Reb_True'.replace('/', ':')
fnames = [one, two, three, four]

plot_comparison(fnames,'4c')




