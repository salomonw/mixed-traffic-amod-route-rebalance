import src.tnet as tnet
import src.CARS as cars
from src.utils import *
import numpy as np
import copy
import matplotlib as mpl
from matplotlib import rc
import matplotlib.pyplot as plt
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)



#netFile, gFile, fcoeffs, tstamp, dir_out = tnet.get_network_parameters('Braess1')
netFile, gFile, fcoeffs, tstamp, dir_out = tnet.get_network_parameters('EMA', experiment_name='EMA_penRate_CARS')
#netFile, gFile, fcoeffs, tstamp, dir_out = tnet.get_network_parameters('NYC_small')
#netFile, gFile, fcoeffs, tstamp, dir_out = tnet.get_network_parameters('NYC_Uber_small', experiment_name='NYC_Uber_small_penRate_comparison')
xa = 0.8




print('Penetration Rate Experiment')
print('------------------------------------------------------------------------------------')
print('PenRate \t AMoD \t Private \t Total')
print('------------------------------------------------------------------------------------')
cavsCost = []
noncavsCost = []
totCost = []
cavsFlow = []
nonCavsFlow = []
pedestrianFlow = []
rebalancingFlow = []
for penetration_rate in np.linspace(0.01,0.99, 11):
    tNet_cavs = tnet.tNet(netFile=netFile, gFile=gFile, fcoeffs=fcoeffs)
    tNet_non_cavs = copy.deepcopy(tNet_cavs)
    g_cavs = tnet.perturbDemandConstant(tNet_cavs.g, constant=penetration_rate/37)
    g_non_cavs = tnet.perturbDemandConstant(tNet_cavs.g, constant=(1-penetration_rate)/37)
    tNet_cavs.set_g(g_cavs)
    tNet_non_cavs.set_g(g_non_cavs)
    tNet_cavs.build_supergraph()
    tNet_non_cavs.build_supergraph()

    it = []
    for i in range(7):
        if i==0:
            #tNet_non_cavs.solveMSA()
            cars.solve_CARS2(tNet_cavs, exogenous_G=False, fcoeffs=fcoeffs, rebalancing=False)
           # cars.solve_CARS(tNet_cavs, exogenous_G=False, fcoeffs=fcoeffs, xa=xa, rebalancing=False)
        else:
            #tNet_non_cavs.solveMSA(exogenous_G=tNet_cavs.G_supergraph)
            cars.solve_CARS2(tNet_cavs, exogenous_G=tNet_non_cavs.G_supergraph, fcoeffs=fcoeffs, rebalancing=False)
            #cars.solve_CARS(tNet_cavs, exogenous_G=tNet_non_cavs.G_supergraph, fcoeffs=fcoeffs, xa=xa, rebalancing=False)

        cars.supergraph2G(tNet_cavs)
        tNet_non_cavs.solveMSA(exogenous_G=tNet_cavs.G_supergraph)
        cars.G2supergraph(tNet_non_cavs)

        totalCost = tnet.get_totalTravelTime(tNet_non_cavs.G, fcoeffs, G_exogenous=tNet_cavs.G) + cars.get_totalTravelTime_without_Rebalancing(tNet_cavs,  G_exogenous=tNet_non_cavs.G)
        it.append(totalCost)

    if penetration_rate>.48 and penetration_rate < 0.52:
        mkdir_n('results/' + dir_out)
        plt.figure(num=None, figsize=(3.8, 5))
        plt.plot(it)
        plt.xlabel('Iteration')
        plt.ylabel('Total Objective')
        plt.tight_layout()
        plt.savefig('results/' + dir_out + '/iteration.png')



    cavCost = cars.get_totalTravelTime_without_Rebalancing(tnet=tNet_cavs, G_exogenous=tNet_non_cavs.G_supergraph)
    nonCavCost = tnet.get_totalTravelTime(tNet_non_cavs.G_supergraph, fcoeffs=fcoeffs, G_exogenous=tNet_cavs.G_supergraph)
    total_Cost = cavCost+nonCavCost

    cavsCost.append(cavCost/tNet_cavs.totalDemand*60)
    noncavsCost.append(nonCavCost/tNet_non_cavs.totalDemand*60)
    totCost.append(total_Cost/(tNet_cavs.totalDemand+tNet_non_cavs.totalDemand)*60)  #TODO: calculate the cost with rebalancing

    print(str(round(penetration_rate, 2)) + '\t' + str(round(cavCost/tNet_cavs.totalDemand*60, 2)) + '\t' + str(
        round(nonCavCost/tNet_non_cavs.totalDemand*60, 2)) + '\t' + str(round(total_Cost/(tNet_cavs.totalDemand+tNet_non_cavs.totalDemand)*60, 2)))

    cavsFlow.append(cars.get_amod_flow(tNet_cavs))
    nonCavsFlow.append(tnet.get_total_G_flow(tNet_non_cavs.G))
    pedestrianFlow.append(cars.get_pedestrian_flow(tNet_cavs))
    rebalancingFlow.append(cars.get_rebalancing_flow(tNet_cavs))

    del tNet_cavs, tNet_non_cavs


mkdir_n('results/' + dir_out)
mpl.rc('font',**{'family':'Times New Roman', 'size': 12})
plt.figure(num=None, figsize=(3.8, 5))
plt.plot(list(np.linspace(0.0001,1, 11)), totCost, label='Total')
plt.plot(list(np.linspace(0.0001,1, 11)), cavsCost, label='AMoD')
plt.plot(list(np.linspace(0.0001,1, 11)), noncavsCost, label='Private Vehicles')
plt.legend()
plt.xlabel('Penetration Rate')
plt.ylabel('Avg. Travel Time (min)')
plt.tight_layout()
plt.savefig('results/' + dir_out +'/costs.png')


plt.figure(num=None, figsize=(5, 2.8))
width = 0.5
ind = np.arange(11)
p1 = plt.bar(ind, nonCavsFlow, width)
p2 = plt.bar(ind, cavsFlow, width,
             bottom=nonCavsFlow)
p3 = plt.bar(ind, rebalancingFlow, width,
             bottom=[cavsFlow[i] + nonCavsFlow[i] for i in range(len(cavsFlow))])
p4 = plt.bar(ind, pedestrianFlow, width,
             bottom=[cavsFlow[i] + nonCavsFlow[i] + rebalancingFlow[i] for i in range(len(cavsFlow))])
plt.ylabel('Flow')
plt.xlabel('Penetration Rate')
plt.legend((p1[0], p2[0], p3[0], p4[0]), ('Private Vehicles', 'AMoDs', 'Rebalancing', 'Pedestrian'))
plt.tight_layout()

#plt.show()


plt.savefig('results/' + dir_out +'/flow_composition.png')
plt.show()


#TODO: check social solution and if they are solving for their own fleet