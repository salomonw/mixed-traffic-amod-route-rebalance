import src.tnet as tnet
import networkx as nx
import src.CARS as cars
import matplotlib.pyplot as plt
import numpy as np


def PoA_experiment(netFile, gFile, posFile, fcoeffs, walk_multiplier, xa):
    priceOfAnarchy = []
    percentagePed = []
    x = []
    exo = []
    amod = []
    amod_flow = []

    for g_multiplier in np.linspace(0.025, 4, 25):
        print(g_multiplier)

        tNet = tnet.tNet(netFile=netFile, gFile=gFile, fcoeffs=fcoeffs)
        tNet.read_node_coordinates(posFile)
        tNet.build_supergraph(walk_multiplier=walk_multiplier)

        pedestrian = [(u, v) for (u, v, d) in tNet.G_supergraph.edges(data=True) if d['type'] == 'p']
        connector = [(u, v) for (u, v, d) in tNet.G_supergraph.edges(data=True) if d['type'] == 'f']

        g_k = tnet.perturbDemandConstant(tNet.g.copy(), constant=g_multiplier)
        tNet.set_g(g_k)
        tNet = cars.solve_CARS(tNet, exogenous_G=0, fcoeffs=fcoeffs, xa=xa, rebalancing=False)
        tNet.solveMSA()
        exogObj = tnet.get_totalTravelTime(tNet.G, fcoeffs)
        amodObjNoRebalancing = cars.get_totalTravelTime(tNet)
        exo.append(exogObj)
        amod.append(amodObjNoRebalancing)
        priceOfAnarchy.append(exogObj / amodObjNoRebalancing)
        percentagePed.append(
            cars.get_pedestrian_flow(tNet) / (cars.get_pedestrian_flow(tNet) + cars.get_amod_flow(tNet)) * 100)
        amod_flow.append(cars.get_amod_flow(tNet))
        x.append(g_multiplier)
    return tNet, priceOfAnarchy, percentagePed, amod_flow, x



netFile, gFile, fcoeffs = tnet.get_network_parameters('Braess1')
posFile = 'data/pos/Braess1_pos.txt'
xa = 0.01

fig, ax1 = plt.subplots()

tNet, priceOfAnarchy, percentagePed, amod_flow, x = PoA_experiment(netFile, gFile, posFile, fcoeffs, walk_multiplier=100000000, xa=xa)

tnet.plot_network_flows(tNet.G, width=3, cmap=plt.cm.Blues)
cars.plot_supergraph_car_flows(tNet)

fig1, ax1 = plt.subplots()
ax1.plot(x, priceOfAnarchy, '--',  label='Price of Anarchy no walking', color='black')
#ax1.legend(loc=1)

tNet, priceOfAnarchy, percentagePed, amod_flow, x = PoA_experiment(netFile, gFile, posFile, fcoeffs, walk_multiplier=4, xa=xa)
ax1.plot(x, priceOfAnarchy, label='Price of Anarchy', color='red')
ax1.legend(loc=0)


ax2 = ax1.twinx()
ax2.plot(x, percentagePed, label='% of Pedestrians', color='blue')
ax2.plot(x, [i/amod_flow[-1:][0]*100 for i in amod_flow], label='Total AMoD veh flow', color='green')
ax2.legend(loc=1)
ax2.set_xlabel('Demand multiplier')
ax1.set_ylabel('PoA')
ax2.set_ylabel('(%)')
#fig1.tight_layout()


tnet.plot_network_flows(tNet.G, width=3, cmap=plt.cm.Blues)

plt.show()