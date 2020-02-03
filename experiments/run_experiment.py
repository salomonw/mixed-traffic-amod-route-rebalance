import sys
from src.utils import *
import src.tnet as tnet
import os
import networkx as nx
import matplotlib.pyplot as plt
import src.CARS as cars

netFile = "data/net/EMA_net.txt"
gFile = "data/trips/EMA_trips.txt"
fcoeffs = [1,0,0,0,0.15,0]

tNet = tnet.tNet(netFile=netFile, gFile=gFile, fcoeffs=fcoeffs)
tNetExog = tnet.tNet(netFile=netFile, gFile=gFile, fcoeffs=fcoeffs)

g_exog = tnet.perturbDemandConstant(tNetExog.g, constant=1.5)
tNetExog.set_g(g_exog)

tNet.build_supergraph()

pedestrian = [(u, v) for (u, v, d) in tNet.G_supergraph.edges(data=True) if d['type'] == 'p']
connector = [(u, v) for (u, v, d) in tNet.G_supergraph.edges(data=True) if d['type'] == 'f']


totalCost = []
for i in range(5):
    if i == 0:
        tNetExog.solveMSA()
    else:
        tNetExog.solveMSA(exogenous_G=tNet.G_supergraph)
    cars.solve_CARS2(tNet, exogenous_G=tNetExog.G, fcoeffs=fcoeffs)
    G_total = tnet.add_G_flows([tNetExog.G, tNet.G_supergraph])
    totalCost.append(tnet.get_totalTravelTime(G_total, fcoeffs))

plt.plot(totalCost)
plt.show()