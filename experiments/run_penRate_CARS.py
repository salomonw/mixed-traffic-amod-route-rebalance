
import src.tnet as tnet
import src.CARS as cars
import numpy as np
import matplotlib.pyplot as plt
import copy

netFile = "data/net/NYC_small_net.txt"
gFile = "data/trips/NYC_small_trips.txt"
fcoeffs = [1,0,0,0,0.15,0]


cavsCost=[]
noncavsCost=[]
totCost = []
for penetration_rate in np.linspace(0.01,0.99, 10):
    tNet_cavs = tnet.tNet(netFile=netFile, gFile=gFile, fcoeffs=fcoeffs)
    tNet_cavs.build_supergraph(walk_multiplier=9999999.12)
    tNet_non_cavs = copy.deepcopy(tNet_cavs)
    g_cavs = tnet.perturbDemandConstant(tNet_cavs.g, constant=penetration_rate*1)
    g_non_cavs = tnet.perturbDemandConstant(tNet_cavs.g, constant=(1-penetration_rate)*1)
    tNet_cavs.set_g(g_cavs)
    tNet_non_cavs.set_g(g_non_cavs)


    it = []
    for i in range(5):
        if i==0:
            tNet_non_cavs.solveMSA()
        else:
            tNet_non_cavs.solveMSA(exogenous_G=tNet_cavs.G_supergraph)

        #cars.solve_CARS2(tNet_cavs, exogenous_G=tNet_non_cavs.G, fcoeffs=fcoeffs)
        cars.solve_CARS2_noRebalancing(tNet_cavs, exogenous_G=tNet_non_cavs.G, fcoeffs=fcoeffs)
        cars.supergraph2G(tNet_cavs)
        totalCost = tnet.get_totalTravelTime(tNet_non_cavs.G,fcoeffs, G_exogenous=tNet_cavs.G) + cars.get_totalTravelTime_without_Rebalancing(tNet_cavs, G_exogenous=tNet_non_cavs.G)
        it.append(totalCost)


    if penetration_rate>.4 and penetration_rate < 0.58:
        plt.figure()
        plt.plot(it)


    cavsCost.append(cars.get_totalTravelTime_without_Rebalancing(tNet_cavs, G_exogenous=tNet_non_cavs.G)/tNet_cavs.totalDemand)
    noncavsCost.append(tnet.get_totalTravelTime(tNet_non_cavs.G,fcoeffs, G_exogenous=tNet_cavs.G)/tNet_non_cavs.totalDemand)
    totCost.append(totalCost/(tNet_cavs.totalDemand+tNet_non_cavs.totalDemand))  #TODO: calculate the cost with rebalancing
    print(penetration_rate)
    del tNet_cavs, tNet_non_cavs


plt.figure()
plt.plot(list(np.linspace(0.0001,1, 10)), totCost, label='Total, t(x+u+r)(x+u)')
plt.plot(list(np.linspace(0.0001,1, 10)), cavsCost, label='CAV, t(x+u+r)(x)')
plt.plot(list(np.linspace(0.0001,1, 10)), noncavsCost, label='NonCAV, t(x+u+r)(u)')
plt.legend()
plt.show()

#TODO: check social solution and if they are solving for their own fleet