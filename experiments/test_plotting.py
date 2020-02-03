import src.tnet as tnet
import networkx as nx
import src.CARS as cars
import matplotlib.pyplot as plt
import numpy as np


netFile = "data/net/Braess1_net.txt"
gFile = "data/trips/Braess1_trips.txt"
posFile = "data/pos/Braess1_pos.txt"

fcoeffs = [1,1,0,0,0,0]
tNet = tnet.tNet(netFile=netFile, gFile=gFile, fcoeffs=fcoeffs)

tNet.read_node_coordinates(posFile)
tNet.solveMSA()

pos = nx.get_node_attributes(tNet.G,'pos')
nx.draw(tNet.G, pos)
plt.show()

tNet.build_supergraph(walk_multiplier=walk_multiplier)
