

import src.tnet as tnet
import numpy as np
import matplotlib.pyplot as plt
import copy

netFile = "data/net/EMA_net.txt"
gFile = "data/trips/EMA_trips.txt"
fcoeffs = [1,0,0,0,0.15,0]


cavsCost=[]
noncavsCost=[]
totCost = []
for penetration_rate in np.linspace(0.01,0.99, 10):
    tNet_cavs = tnet.tNet(netFile=netFile, gFile=gFile, fcoeffs=fcoeffs)
    tNet_non_cavs = copy.deepcopy(tNet_cavs)
    g_cavs = tnet.perturbDemandConstant(tNet_cavs.g, constant=penetration_rate)
    g_non_cavs = tnet.perturbDemandConstant(tNet_cavs.g, constant=(1-penetration_rate))
    tNet_cavs.set_g(g_cavs)
    tNet_non_cavs.set_g(g_non_cavs)


    it = []
    for i in range(5):
        if i==0:
            # tnet.solveMSA_julia_social(tNet_cavs)
            tNet_cavs.solveMSAsocial()
        else:
            # tnet.solveMSA_julia_social(tNet_cavs, G_exo=tNet_non_cavs.G)
            tNet_cavs.solveMSAsocial(exogenous_G=tNet_non_cavs.G)

        #tnet.solveMSA_julia(tNet_non_cavs, G_exo=tNet_cavs.G)
        tNet_non_cavs.solveMSA(exogenous_G=tNet_cavs.G)
        G_total = tnet.add_net_flows([tNet_non_cavs, tNet_cavs])
        totalCost = tnet.get_totalTravelTime(G_total, fcoeffs)
        it.append(totalCost)

    if penetration_rate>.48 and penetration_rate < 0.58:
        plt.figure()
        plt.plot(it)

    cavsCost.append(tnet.get_totalTravelTime(tNet_cavs.G, fcoeffs, G_exogenous=tNet_non_cavs.G)/tNet_cavs.totalDemand*60)
    noncavsCost.append(tnet.get_totalTravelTime(tNet_non_cavs.G,fcoeffs, G_exogenous=tNet_cavs.G)/tNet_non_cavs.totalDemand*60)

    G_total = tnet.add_net_flows([tNet_non_cavs, tNet_cavs])
    totalCost = tnet.get_totalTravelTime(G_total, fcoeffs)
    totCost.append(tnet.get_totalTravelTime(G_total, fcoeffs)/(tNet_cavs.totalDemand+tNet_non_cavs.totalDemand)*60)
    print(penetration_rate)
    del tNet_cavs, tNet_non_cavs

plt.figure()
plt.plot(list(np.linspace(0.0001,1, 10)), totCost, label='Total')
plt.plot(list(np.linspace(0.0001,1, 10)), cavsCost, label='CAV')
plt.plot(list(np.linspace(0.0001,1, 10)), noncavsCost, label='NonCAV')
plt.legend()
plt.show()

#TODO: check social solution and if they are solving for their own fleet