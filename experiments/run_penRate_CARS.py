import src.tnet as tnet
import src.CARS as cars
import numpy as np
import matplotlib.pyplot as plt
import copy


#netFile, gFile, fcoeffs = tnet.get_network_parameters('Braess1')
netFile, gFile, fcoeffs = tnet.get_network_parameters('EMA')
#netFile, gFile, fcoeffs = tnet.get_network_parameters('NYC_small')
xa = 0.8

cavsCost = []
noncavsCost = []
totCost = []
cavsFlow = []
nonCavsFlow = []
pedestrianFlow = []
rebalancingFlow = []
for penetration_rate in np.linspace(0.01,0.99, 10):
    tNet_cavs = tnet.tNet(netFile=netFile, gFile=gFile, fcoeffs=fcoeffs)
    tNet_non_cavs = copy.deepcopy(tNet_cavs)
    g_cavs = tnet.perturbDemandConstant(tNet_cavs.g, constant=penetration_rate)
    g_non_cavs = tnet.perturbDemandConstant(tNet_cavs.g, constant=(1-penetration_rate))
    tNet_cavs.set_g(g_cavs)
    tNet_non_cavs.set_g(g_non_cavs)
    tNet_cavs.build_supergraph(walk_multiplier=.80)

    it = []
    for i in range(10):
        if i==0:
            #tNet_non_cavs.solveMSA()
            #cars.solve_CARS2(tNet_cavs, exogenous_G=False, fcoeffs=fcoeffs, xa=xa, rebalancing=False)
            cars.solve_CARS(tNet_cavs, exogenous_G=False, fcoeffs=fcoeffs, xa=xa, rebalancing=True)
        else:
            #tNet_non_cavs.solveMSA(exogenous_G=tNet_cavs.G_supergraph)
            #cars.solve_CARS2(tNet_cavs, exogenous_G=tNet_non_cavs.G, fcoeffs=fcoeffs, xa=xa, rebalancing=False)
            cars.solve_CARS(tNet_cavs, exogenous_G=tNet_non_cavs.G, fcoeffs=fcoeffs, xa=xa, rebalancing=True)

        #cars.solve_CARS2(tNet_cavs, exogenous_G=tNet_non_cavs.G, fcoeffs=fcoeffs, xa=xa, rebalancing=False)
        #cars.solve_CARS(tNet_cavs, exogenous_G=tNet_non_cavs.G, fcoeffs=fcoeffs, xa=xa, rebalancing=False)
        cars.supergraph2G(tNet_cavs)
        tNet_non_cavs.solveMSA(exogenous_G=tNet_cavs.G)

        totalCost = tnet.get_totalTravelTime(tNet_non_cavs.G, fcoeffs, G_exogenous=tNet_cavs.G) + cars.get_totalTravelTime_without_Rebalancing(tNet_cavs,  G_exogenous=tNet_non_cavs.G)
        it.append(totalCost)


    if penetration_rate>.4 and penetration_rate < 0.58:
        plt.figure()
        plt.plot(it)


    cavsCost.append(cars.get_totalTravelTime_without_Rebalancing(tNet_cavs, G_exogenous=tNet_non_cavs.G)/tNet_cavs.totalDemand)
    noncavsCost.append(tnet.get_totalTravelTime(tNet_non_cavs.G,fcoeffs, G_exogenous=tNet_cavs.G)/tNet_non_cavs.totalDemand)
    totCost.append(totalCost/(tNet_cavs.totalDemand+tNet_non_cavs.totalDemand))  #TODO: calculate the cost with rebalancing

    cavsFlow.append(cars.get_amod_flow(tNet_cavs))
    nonCavsFlow.append(tnet.get_total_G_flow(tNet_non_cavs.G))
    pedestrianFlow.append(cars.get_pedestrian_flow(tNet_cavs))
    rebalancingFlow.append(cars.get_rebalancing_flow(tNet_cavs))

    print(penetration_rate)
    del tNet_cavs, tNet_non_cavs


plt.figure()
plt.plot(list(np.linspace(0.0001,1, 10)), totCost, label='Total, t(x+u+r)(x+u)')
plt.plot(list(np.linspace(0.0001,1, 10)), cavsCost, label='CAV, t(x+u+r)(x)')
plt.plot(list(np.linspace(0.0001,1, 10)), noncavsCost, label='NonCAV, t(x+u+r)(u)')
plt.legend()

plt.figure()
width = 0.5
ind = np.arange(10)
p1 = plt.bar(ind, nonCavsFlow, width)
p2 = plt.bar(ind, cavsFlow, width,
             bottom=nonCavsFlow)
p3 = plt.bar(ind, rebalancingFlow, width,
             bottom=[cavsFlow[i] + nonCavsFlow[i] for i in range(len(cavsFlow))])
p4 = plt.bar(ind, pedestrianFlow, width,
             bottom=[cavsFlow[i] + nonCavsFlow[i] + rebalancingFlow[i] for i in range(len(cavsFlow))])
plt.ylabel('Flow')
plt.legend((p1[0], p2[0], p3[0], p4[0]), ('NonCAVs', 'CAVs', 'Rebalancing', 'Pedestrian'))
plt.show()



#TODO: check social solution and if they are solving for their own fleet