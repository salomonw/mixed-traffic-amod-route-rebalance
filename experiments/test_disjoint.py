from src.utils import *
import src.tnet as tnet
import matplotlib.pyplot as plt
import src.CARS as cars

#netFile, gFile, fcoeffs = tnet.get_network_parameters('Braess1')
netFile, gFile, fcoeffs = tnet.get_network_parameters('EMA')
#netFile, gFile, fcoeffs = tnet.get_network_parameters('NYC_small')

tNet = tnet.tNet(netFile=netFile, gFile=gFile, fcoeffs=fcoeffs)
tNet.build_supergraph()

pedestrian = [(u, v) for (u, v, d) in tNet.G_supergraph.edges(data=True) if d['type'] == 'p']
connector = [(u, v) for (u, v, d) in tNet.G_supergraph.edges(data=True) if d['type'] == 'f']

#TODO: MSA social be able to solve with pedestrian network, implement function in CARS modusle calling tnet function
tNet.solveMSAsocial()
cars.solve_rebalancing(tNet, exogenous_G=0)

