import src.tnet as tnet
import src.CARS as cars
import numpy as np
import matplotlib.pyplot as plt
import copy
from src.utils import *

#netFile, gFile, fcoeffs = tnet.get_network_parameters('Braess1')
netFile, gFile, fcoeffs, tstamp, dir_out = tnet.get_network_parameters('EMA', experiment_name='EMA_test_ITSC_2019')
netFile, gFile, fcoeffs, tstamp, dir_out = tnet.get_network_parameters('NYC_Uber_small', experiment_name='NYC_Uber_test_ITSC_2019')
#netFile, gFile, fcoeffs, tstamp, dir_out = tnet.get_network_parameters('Braess1', experiment_name='Braess1_test_ITSC_2019')
#netFile, gFile, fcoeffs = tnet.get_network_parameters('NYC_small')
#netFile, gFile, fcoeffs = tnet.get_network_parameters('NYC_Uber_small')

print('---------------------------------------------')
print('Pen Rate \t POA')
print('---------------------------------------------')
cavsCost=[]
noncavsCost=[]
totCost = []
for penetration_rate in np.linspace(0.01,0.99, 11):
    tNet_cavs = tnet.tNet(netFile=netFile, gFile=gFile, fcoeffs=fcoeffs)
    tNet_cavs.build_supergraph()
    pedestrian = [(u, v) for (u, v, d) in tNet_cavs.G_supergraph.edges(data=True) if d['type'] == 'p']
    connector = [(u, v) for (u, v, d) in tNet_cavs.G_supergraph.edges(data=True) if d['type'] == 'f']

    tNet_non_cavs = copy.deepcopy(tNet_cavs)
    g_cavs = tnet.perturbDemandConstant(tNet_cavs.g, constant=penetration_rate/38)
    g_non_cavs = tnet.perturbDemandConstant(tNet_cavs.g, constant=(1-penetration_rate)/38)
    tNet_cavs.set_g(g_cavs)
    tNet_non_cavs.set_g(g_non_cavs)
    cars.G2supergraph(tNet_non_cavs)
    cars.G2supergraph(tNet_cavs)

    it = []
    for i in range(6):
        if i==0:
            #tNet_cavs.solveMSAsocial()
            cars.solveMSAsocialCARS(tNet_cavs, exogenous_G=False)
            #cars.solve_social_Julia(tNet_cavs, exogenous_G=False)

        else:
            #tnet.solveMSA_julia_social(tNet_cavs, G_exo=tNet_non_cavs.G)
            #cars.solve_social_Julia(tNet_cavs, exogenous_G=tNet_non_cavs.G_supergraph)
            cars.solveMSAsocialCARS(tNet_cavs, exogenous_G=tNet_non_cavs.G_supergraph)

        tNet_non_cavs.solveMSA(exogenous_G=tNet_cavs.G_supergraph)

    if penetration_rate>.48 and penetration_rate < 0.58:
        plt.figure()
        plt.plot(it)
        mkdir_n('results/' + dir_out)
        plt.savefig('results/' + dir_out + '/iteration.png', dpi=300)


    cars.G2supergraph(tNet_non_cavs)

    cavs_cost = tnet.get_totalTravelTime(tNet_cavs.G_supergraph, fcoeffs,
                                         G_exogenous=tNet_non_cavs.G_supergraph)
    noncavs_cost = tnet.get_totalTravelTime(tNet_non_cavs.G_supergraph, fcoeffs, G_exogenous=tNet_cavs.G_supergraph)
    totalCost = cavs_cost + noncavs_cost

    cavsCost.append(cavs_cost/tNet_cavs.totalDemand*60)
    noncavsCost.append(noncavs_cost/tNet_non_cavs.totalDemand*60)
    totCost.append(totalCost/(tNet_cavs.totalDemand+tNet_non_cavs.totalDemand)*60)

    print(str(round(penetration_rate,2)) + '\t' + str(round(totCost[0]/totCost[-1],4)))

    del tNet_cavs, tNet_non_cavs

mkdir_n('results/' + dir_out)
plt.savefig('results/' + dir_out +'/iteration.png', dpi=300)

plt.figure()
plt.plot(list(np.linspace(0.0001,1, 11)), totCost, label='Total')
plt.plot(list(np.linspace(0.0001,1, 11)), cavsCost, label='AMoD')
plt.plot(list(np.linspace(0.0001,1, 11)), noncavsCost, label='Private Veh')
plt.legend()
#plt.show()

plt.savefig('results/' + dir_out +'/penetration_rate.png', dpi=300)
