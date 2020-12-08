import src.tnet as tnet
import src.CARS as cars
import src.routeFinder as routeFinder
import experiments.build_NYC_subway_net as nyc



#netFile, gFile, fcoeffs, tstamp, dir_out = tnet.get_network_parameters('EMA', experiment_name='EMA_penRate_comparison-'+'REB')
#netFile, gFile, fcoeffs, tstamp, dir_out = tnet.get_network_parameters('NYC_Uber_small_1', experiment_name='NYC_Uber_small_1_penRate_comparison-REB')
#netFile, gFile, fcoeffs, tstamp, dir_out = tnet.get_network_parameters('Anaheim', experiment_name='Anaheim_test_CARSn')
#netFile, gFile, fcoeffs, tstamp, dir_out = tnet.get_network_parameters('Braess1', experiment_name='Braess1_test_rebalancingRoutes')
#netFile, gFile, fcoeffs, tstamp, dir_out = tnet.get_network_parameters('Sioux Falls', experiment_name='Sioux-Falls_test_rebalancingRoutes')
#netFile, gFile, fcoeffs, tstamp, dir_out = tnet.get_network_parameters('EMA_mid', experiment_name='EMA_mid_test_rebalancingRoutes')
netFile, gFile, fcoeffs, tstamp, dir_out = tnet.get_network_parameters('EMA_mid', experiment_name='EMA_mid_test_rebalancingRoutes')
#tNet, tstamp, fcoeffs = nyc.build_NYC_net('data/net/NYC_M/', only_road=True)

tNet = tnet.tNet(netFile=netFile, gFile=gFile, fcoeffs=fcoeffs)
tNet.build_supergraph(walk_multiplier=1)

tNet, runtime, od_flows = cars.solve_bush_CARSn(tnet=tNet, fcoeffs=fcoeffs,
                                                n=6, rebalancing=True, linear=False,
                                                bush=True)

print('CARS solved!')
cars.supergraph2G(tNet)
tNet.set_supergraph_tk(fcoeffs=fcoeffs, G_exogenous=False)
rebRoutes = routeFinder.rebRouteFinder(tNet.G, eps=10, print_=True)
print('Rebalanacing routes found!')
print(rebRoutes)

