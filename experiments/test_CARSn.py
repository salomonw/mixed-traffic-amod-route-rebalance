import src.tnet as tnet
import src.CARS as cars
import numpy as np
import copy
from src.utils import *
import matplotlib as mpl
from matplotlib import rc
import matplotlib.pyplot as plt
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#rc('text', usetex=True)

netFile, gFile, fcoeffs, tstamp, dir_out = tnet.get_network_parameters('EMA', experiment_name='EMA_penRate_comparison-'+'REB')
#netFile, gFile, fcoeffs, tstamp, dir_out = tnet.get_network_parameters('NYC_Uber_small', experiment_name='NYC_Uber_small_penRate_comparison')
#netFile, gFile, fcoeffs, tstamp, dir_out = tnet.get_network_parameters('NYC_Uber_small_1', experiment_name='NYC_Uber_small_1_penRate_comparison-REB')
#netFile, gFile, fcoeffs, tstamp, dir_out = tnet.get_network_parameters('Anaheim', experiment_name='Anaheim_test_CARSn')
#netFile, gFile, flowFile, fcoeffs, tstamp, dir_out = tnet.get_network_parameters('Barcelons', experiment_name='Barcelona_buildNet')
#netFile, gFile, flowFile, fcoeffs, tstamp, dir_out = tnet.get_network_parameters('ChicagoSketch', experiment_name='ChicagoSketch')
#netFile, gFile, flowFile, fcoeffs, tstamp, dir_out = tnet.get_network_parameters('Sydeny', experiment_name='Sydeny')

demand_multiplier = list(np.linspace(0.8,1.8,2))
demand_multiplier = [1]

'''

print('---- solving NLP problem to set up a base ---')
real_obj = []
for g_multi in demand_multiplier:
    tNet = tnet.tNet(netFile=netFile, gFile=gFile, fcoeffs=fcoeffs)
    tNet.build_supergraph(walk_multiplier=1)
    pedestrian = [(u, v) for (u, v, d) in tNet.G_supergraph.edges(data=True) if d['type'] == 'p']
    connector = [(u, v) for (u, v, d) in tNet.G_supergraph.edges(data=True) if d['type'] == 'f']
    g_per = tnet.perturbDemandConstant(tNet.g, g_multi)
    tNet.set_g(g_per)
    cars.solve_social_Julia(tNet, exogenous_G=False)
    print('\t solve for g_multiplier = ' + str(round(g_multi,2)))
    socialObj = tnet.get_totalTravelTime(tNet.G_supergraph, fcoeffs)
    real_obj.append(socialObj)
    print(socialObj)
'''


n = [2+i for i in range(4)]

print("\ntestCars progressBar:")
progBar = progressBar(len(n)*2*len(demand_multiplier))
progBar.set()
CARS = {}
for i in n:
    CARS[i] = {}
    for g_multi in demand_multiplier:
        for linear in [True, False]:
            tNet = tnet.tNet(netFile=netFile, gFile=gFile, fcoeffs=fcoeffs)
            tNet.build_supergraph(walk_multiplier=1)
            pedestrian = [(u, v) for (u, v, d) in tNet.G_supergraph.edges(data=True) if d['type'] == 'p']
            connector = [(u, v) for (u, v, d) in tNet.G_supergraph.edges(data=True) if d['type'] == 'f']

            g_per = tnet.perturbDemandConstant(tNet.g, g_multi)
            tNet.set_g(g_per)

            tNet, runtime, od_flows = cars.solve_CARSn(tNet, fcoeffs=fcoeffs, n=i, exogenous_G=False, rebalancing=False, linear=linear,  method=1)

            CARS2obj = tnet.get_totalTravelTime(tNet.G_supergraph, fcoeffs)
            CARS[i][linear] = (CARS2obj-1630.1380990494615)/1630.1380990494615*100

            progBar.tic()
        del tNet

fig, ax = plt.subplots(figsize=(5,2))
ax.plot(n, [v[True] for k,v in CARS.items()], label = 'LP')
ax.plot(n, [v[False] for k,v in CARS.items()], label = 'QP')
ax.set_xlabel('n')
ax.set_ylabel('% deviation from NLP')
ax.set_xlim([n[0], n[-1]])
ax.legend(framealpha=1)
ax.grid(True)
#plt.tight_layout()
plt.show()


