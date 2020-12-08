from src.utils import *
import src.tnet as tnet
import src.CARS as cars
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

#netFile, gFile, fcoeffs, tstamp, dir_out = tnet.get_network_parameters('EMA', experiment_name='EMA_CARS_precision')
#netFile, gFile, fcoeffs, tstamp, dir_out = tnet.get_network_parameters('NYC_Uber_small_1', experiment_name='NYC_Uber_small_1_precision')
netFile, gFile, fcoeffs, tstamp, dir_out = tnet.get_network_parameters('Anaheim', experiment_name='Anaheim_poa_finder')

demand_multiplier = list(np.linspace(0.5,2.5,11))

for g_multi in demand_multiplier:
    # BUILD NETWORK
    tNet = tnet.tNet(netFile=netFile, gFile=gFile, fcoeffs=fcoeffs)
    tNet.build_supergraph()
    pedestrian = [(u, v) for (u, v, d) in tNet.G_supergraph.edges(data=True) if d['type'] == 'p']
    connector = [(u, v) for (u, v, d) in tNet.G_supergraph.edges(data=True) if d['type'] == 'f']
    g_per = tnet.perturbDemandConstant(tNet.g, g_multi)
    tNet.set_g(g_per)

    #if g_multi>2.65:
    #    tnet.writeTripsfile(tNet.g, 'hola2.txt')
    #    dsd
    # SYSTEM-CENTRIC
    tNet.solveMSAsocial_supergraph()
    cars.supergraph2G(tNet)
    sc = cars.get_totalTravelTime(tNet)

    del tNet
    # BUILD NETWORK
    tNet = tnet.tNet(netFile=netFile, gFile=gFile, fcoeffs=fcoeffs)
    tNet.build_supergraph()
    pedestrian = [(u, v) for (u, v, d) in tNet.G_supergraph.edges(data=True) if d['type'] == 'p']
    connector = [(u, v) for (u, v, d) in tNet.G_supergraph.edges(data=True) if d['type'] == 'f']
    g_per = tnet.perturbDemandConstant(tNet.g, g_multi)
    tNet.set_g(g_per)

    # USER-CENTRIC
    tNet.solveMSA()
    cars.G2supergraph(tNet)
    uc = cars.get_totalTravelTime(tNet)

    # PRINT POA
    print('\t g_multi='+ str(round(g_multi,1)) + '\t POA = ' + str(round(uc/sc,4)))