import src.tnet as tnet
import src.CARS as cars
import src.routeFinder as routeFinder

#netFile, gFile, fcoeffs, tstamp, dir_out = tnet.get_network_parameters('EMA', experiment_name='EMA_penRate_comparison-'+'REB')
#netFile, gFile, fcoeffs, tstamp, dir_out = tnet.get_network_parameters('NYC_Uber_small_1', experiment_name='NYC_Uber_small_1_penRate_comparison-REB')
#netFile, gFile, fcoeffs, tstamp, dir_out = tnet.get_network_parameters('Anaheim', experiment_name='Anaheim_test_CARSn')
#netFile, gFile, fcoeffs, tstamp, dir_out = tnet.get_network_parameters('Braess1', experiment_name='Braess1_test_rebalancingRoutes')
#netFile, gFile, fcoeffs, tstamp, dir_out = tnet.get_network_parameters('Sioux Falls', experiment_name='Sioux-Falls_test_rebalancingRoutes')

netFile, gFile, fcoeffs, tstamp, dir_out = tnet.get_network_parameters('EMA_mid', experiment_name='EMA_mid_test_rebalancingRoutes')


tNet = tnet.tNet(netFile=netFile, gFile=gFile, fcoeffs=fcoeffs)
tNet.build_supergraph(walk_multiplier=1)
q = [(u, v) for (u, v, d) in tNet.G_supergraph.edges(data=True) if d['type'] == 'p']
q = [(u, v) for (u, v, d) in tNet.G_supergraph.edges(data=True) if d['type'] == 'f']
g_per = tnet.perturbDemandConstant(tNet.g, 1)
tNet.set_g(g_per)

tNet, runtime, s_flows = cars.solve_bush_CARSn(tNet, fcoeffs=fcoeffs, n=6,
                                                exogenous_G=False, rebalancing=True,
                                                linear=True,  method=1, bush=True)
cars.supergraph2G(tNet)
rebRoutes = routeFinder.rebRouteFinder(tNet.G, eps=0)
userRoutes = routeFinder.userRouteFinder(tNet.G, tNet.g, s_flows, eps=0)
print(userRoutes)


