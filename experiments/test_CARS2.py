import src.tnet as tnet
import src.CARS as cars


netFile, gFile, fcoeffs = tnet.get_network_parameters('NYC_small')


tNet = tnet.tNet(netFile=netFile, gFile=gFile, fcoeffs=fcoeffs)
tNet.build_supergraph(walk_multiplier=0.115)

pedestrian = [(u, v) for (u, v, d) in tNet.G_supergraph.edges(data=True) if d['type'] == 'p']
connector = [(u, v) for (u, v, d) in tNet.G_supergraph.edges(data=True) if d['type'] == 'f']

tNet = cars.solve_CARS2(tNet, exogenous_G=0, fcoeffs=fcoeffs)