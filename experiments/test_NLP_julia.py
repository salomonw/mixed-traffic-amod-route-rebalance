import src.tnet as tnet
import src.CARS as cars
import numpy as np
import matplotlib.pyplot as plt
import copy

#netFile, gFile, fcoeffs = tnet.get_network_parameters('Braess1')
#netFile, gFile, fcoeffs = tnet.get_network_parameters('EMA')
netFile, gFile, fcoeffs = tnet.get_network_parameters('NYC_small')

tNet = tnet.tNet(netFile=netFile, gFile=gFile, fcoeffs=fcoeffs)
tNet.build_supergraph(walk_multiplier=0.115)

pedestrian = [(u, v) for (u, v, d) in tNet.G_supergraph.edges(data=True) if d['type'] == 'p']
connector = [(u, v) for (u, v, d) in tNet.G_supergraph.edges(data=True) if d['type'] == 'f']

print(tNet.G_supergraph.edges(data=True))
tNet.solveMSA()
cars.solve_social_Juliajson(tNet, exogenous_G=tNet.G)

#cars.solve_social_NLP(tNet)