import src.tnet as tnet
import src.CARS as cars
import matplotlib.pyplot as plt
import  numpy as np

netFile, gFile, fcoeffs = tnet.get_network_parameters('EMA')
xa = 0.6

totalObj = []
percentagePed = []
x= []
poa =[]
for g_multiplier in np.linspace(0.025, 3, 20):
    tNet = tnet.tNet(netFile=netFile, gFile=gFile, fcoeffs=fcoeffs)
    tNet.build_supergraph(walk_multiplier=0.115)

    pedestrian = [(u, v) for (u, v, d) in tNet.G_supergraph.edges(data=True) if d['type'] == 'p']
    connector = [(u, v) for (u, v, d) in tNet.G_supergraph.edges(data=True) if d['type'] == 'f']
    g_k = tnet.perturbDemandConstant(tNet.g, constant=g_multiplier)
    totalObj = []
    for j in range(10):
        tNet.set_g(g_k)
        if j==0:
            tNet.solveMSA()
        else:
            tNet.solveMSA(exogenous_G=tNet.G_supergraph)
        tNet = cars.solve_CARS(tNet, exogenous_G=tNet.G, fcoeffs=fcoeffs, xa=xa)
        exogObj = tnet.get_totalTravelTime(tNet.G, fcoeffs)
        amodObj = cars.get_totalTravelTime(tNet)
        amodObjNoRebalancing = cars.get_totalTravelTime_without_Rebalancing(tNet)
        totalObj.append(exogObj + amodObj)
    #plt.plot(totalObj)
    #plt.show()
    print(exogObj/amodObjNoRebalancing)
    poa.append(exogObj/amodObjNoRebalancing)
    percentagePed.append(cars.get_pedestrian_flow(tNet) / (cars.get_pedestrian_flow(tNet)+ cars.get_amod_flow(tNet))*100)
    x.append(g_multiplier)


fig, ax1 = plt.subplots()
ax1.plot(x, poa, '--',  label='Price of Anarchy', color='red')
plt.legend()
ax2 = ax1.twinx()
ax2.plot(x, percentagePed, label='% of Pedestrians', color='blue')
plt.legend()
plt.xlabel('Demand')
plt.ylabel('PoA')
fig.tight_layout()
plt.show()
