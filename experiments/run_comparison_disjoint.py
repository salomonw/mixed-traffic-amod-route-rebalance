import src.tnet as tnet
import src.CARS as cars
from matplotlib import rc
import matplotlib.pyplot as plt
import pandas as pd
import experiments.build_NYC_subway_net as nyc

rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})

def read_net(net_name, g_multi=1):
    if net_name == 'NYC':
        tNet, tstamp, fcoeffs = nyc.build_NYC_net('data/net/NYC_M/', only_road=True)
    else:
        netFile, gFile, fcoeffs, tstamp, dir_out = tnet.get_network_parameters(net_name=net_name,
                                                                               experiment_name=net_name + '_n_variation')
        tNet = tnet.tNet(netFile=netFile, gFile=gFile, fcoeffs=fcoeffs)

    g_per = tnet.perturbDemandConstant(tNet.g, g_multi)
    tNet.set_g(g_per)

    return tNet


def run_comparison_disjoint(net_name, ns=[7], theta_n=3, g_multi=1):
    d = []
    for n in ns:
        for linear in [False]:
            tNet = read_net(net_name, g_multi=g_multi)
            tNet.build_supergraph(walk_multiplier=1, identical_G=True)
            tNet, runtimeCARS, od_flows = cars.solve_bush_CARSn(tNet, fcoeffs=tNet.fcoeffs, n=n, theta_n=theta_n,
                                                            exogenous_G=False, rebalancing=True, linear=linear, bush=True)
            CARSnobj = cars.get_CARS_obj_val(tNet, G_exogenous=False)#(tNet.G_supergraph, tNet.fcoeffs)
            if linear == True:
                tp = 'LP'
            else:
                tp = 'QP'
            d0 = {'Net': net_name, 'A': tNet.nLinks, 'V': tNet.nNodes, 'W': tNet.nOD,
                  'type': 'CARS' + str(n) + '-' + tp, 'obj': CARSnobj, 't': runtimeCARS, 't_reb':0}
            d.append(d0)
            del tNet

    tNet = read_net(net_name, g_multi=g_multi)
    tNet.build_supergraph(walk_multiplier=1, identical_G=True)
    tNet, runtime, od_flows = cars.solve_bush_CARSn(tNet, fcoeffs=tNet.fcoeffs, n=n, theta_n=theta_n,
                                                    exogenous_G=False, rebalancing=False, bush=True, linear=linear,
                                                    userCentric=False)
    #tNet.solveMSAsocial_supergraph()
    cars.supergraph2G(tNet)
    #runtime, RG = cars.solveMSAsocialCARS(tNet, exogenous_G=False)
    runtimeReb = cars.solve_rebalancing(tNet,exogenous_G=0)
    NLP_obj = cars.get_CARS_obj_val(tNet, G_exogenous=False)#tnet.get_totalTravelTime(tNet.G_supergraph, fcoeffs=tNet.fcoeffs, G_exogenous=False)
    d0 = {'Net': net_name, 'A': tNet.nLinks, 'V': tNet.nNodes, 'W': tNet.nOD, 'type': 'dijoint', 'obj': str(round((NLP_obj/CARSnobj-1)*100,2))+"%",
          't': str(round(((runtime+runtimeReb)/runtimeCARS-1)*100,2))+"%", 't_reb':runtimeReb}
    d.append(d0)
    del tNet



    return pd.DataFrame(d)


#nets = ['NYC_Uber_small_1', 'NYC', 'EMA_mid', 'Anaheim', 'Winnipeg', 'Barcelona',  'Sydeny' ]
nets = ["EMA_mid", "NYC"]
#nets = ["EMA_mid"]
#nets = ["EMA"]
#nets = ['Sioux Falls', 'Anaheim', "EMA_mid", "NYC"]
#nets = ['NYC']

ns = [5]

g_multipliers = [2,3]# [.5]#,1,1.5]#,2,4,6]
for net in nets:
    for g in g_multipliers:
        df = run_comparison_disjoint(net, ns=ns, theta_n=g, g_multi=g)
        print(df)
        df.to_csv('comparison_disjoint_scan_' + net + '_' + str(g)+ '.csv')
