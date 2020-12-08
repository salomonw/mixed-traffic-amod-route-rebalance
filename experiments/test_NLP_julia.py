import src.tnet as tnet
import src.CARS as cars
import numpy as np
import matplotlib.pyplot as plt
import copy

#netFile, gFile, fcoeffs = tnet.get_network_parameters('Braess1')
netFile, gFile, fcoeffs, tstamp, dir_out = tnet.get_network_parameters('EMA', experiment_name='EMA_penRate_disjoint')
#netFile, gFile, fcoeffs = tnet.get_network_parameters('NYC_small')

tNet = tnet.tNet(netFile=netFile, gFile=gFile, fcoeffs=fcoeffs)
tNet.build_supergraph()

pedestrian = [(u, v) for (u, v, d) in tNet.G_supergraph.edges(data=True) if d['type'] == 'p']
connector = [(u, v) for (u, v, d) in tNet.G_supergraph.edges(data=True) if d['type'] == 'f']


tNet.solveMSA()
cars.solve_social_Julia(tNet, exogenous_G=False)
#cars.solve_social_Julia(tNet, exogenous_G=tNet.G)

#cars.solve_social_NLP(tNet)