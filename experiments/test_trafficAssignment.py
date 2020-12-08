import src.tnet as tnet
import src.CARS as cars

netFile, gFile, fcoeffs = tnet.get_network_parameters('Braess1')

tNet = tnet.tNet(netFile=netFile, gFile=gFile, fcoeffs=fcoeffs)

tNet.solveMSA()

print([(i,j, tNet.G[i][j]['flow']) for i,j in tNet.G.edges()])

tNet.build_supergraph()
tNet = cars.solve_CARS_noRebalancing(tNet, exogenous_G=0, fcoeffs=fcoeffs, xa=1)

print([(i,j, tNet.G_supergraph[i][j]['flow']) for i,j in tNet.G.edges()])

exogObj = tnet.get_totalTravelTime(tNet.G, fcoeffs)
amodObjNoRebalancing = cars.get_totalTravelTime(tNet)
priceOfAnarchy = exogObj / amodObjNoRebalancing

print(priceOfAnarchy)