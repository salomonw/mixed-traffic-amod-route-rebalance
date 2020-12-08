import src.tnet as tnet
import src.CARS as cars
import os
from datetime import datetime
import experiments.build_NYC_subway_net as nyc


tNet, tstamp, fcoeffs = nyc.build_NYC_net('data/net/NYC_M/')

tNet, runtime, od_flows = cars.solve_bush_CARSn(tNet,
                                           fcoeffs=fcoeffs,
                                           n=4,
                                           exogenous_G=False,
                                           rebalancing=False,
                                           linear=True,
                                           bush=True)

print(tnet.get_totalTravelTime(tNet.G_supergraph, fcoeffs))

#tNet.solveMSAsocial_supergraph(build_t0=False, exogenous_G=False)

print(tnet.get_totalTravelTime(tNet.G_supergraph, fcoeffs))