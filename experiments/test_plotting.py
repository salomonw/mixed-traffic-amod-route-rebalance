import src.tnet as tnet
import networkx as nx
import matplotlib.pyplot as plt

netFile, gFile, fcoeffs = tnet.get_network_parameters('Braess1')
posFile = "data/pos/Braess1_pos.txt"
tNet = tnet.tNet(netFile=netFile, gFile=gFile, fcoeffs=fcoeffs)

tNet.read_node_coordinates(posFile)
tNet.solveMSA()

pos = nx.get_node_attributes(tNet.G,'pos')
nx.draw(tNet.G, pos)
plt.show()

tNet.build_supergraph(walk_multiplier=walk_multiplier)
