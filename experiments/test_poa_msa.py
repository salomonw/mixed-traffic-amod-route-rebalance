import src.tnet as tnet
import numpy as np
import matplotlib.pyplot as plt

netFile, gFile, fcoeffs = tnet.get_network_parameters('EMA')

demand_multiplier = list(np.linspace(0.5, 4, 25))
poa = []
for g_multiplier in demand_multiplier:
    tNet = tnet.tNet(netFile=netFile, gFile=gFile, fcoeffs=fcoeffs)
    g_k = tnet.perturbDemandConstant(tNet.g.copy(), constant=g_multiplier)
    tNet.set_g(g_k)
    #tnet.solveMSA_julia(tNet)
    tNet.solveMSA()
    userCentricCost = tnet.get_totalTravelTime(tNet.G, fcoeffs)

    tNet = tnet.tNet(netFile=netFile, gFile=gFile, fcoeffs=fcoeffs)
    g_k = tnet.perturbDemandConstant(tNet.g.copy(), constant=g_multiplier)
    tNet.set_g(g_k)
    #tnet.solveMSA_julia_social(tNet.G)
    tNet.solveMSAsocial()
    socialCost = tnet.get_totalTravelTime(tNet.G, fcoeffs)


    priceOfAnarchy = userCentricCost / socialCost
    print(priceOfAnarchy)
    poa.append(priceOfAnarchy)

plt.plot(demand_multiplier , poa)
plt.show()