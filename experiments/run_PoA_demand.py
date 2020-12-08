
import src.tnet as tnet
import src.CARS as cars
from matplotlib import rc
import matplotlib.pyplot as plt
import experiments.build_NYC_subway_net as nyc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})


#net_name = 'NYC_Uber_small_1'
#net_name = 'EMA_mid'
net_name = 'NYC'
#net_name = 'Sioux Falls'
#net_name = 'Anaheim'



def read_net(net_name):
    if net_name == 'NYC':
        tNet, tstamp, fcoeffs = nyc.build_NYC_net('data/net/NYC_M/', only_road=True)
    else:
        netFile, gFile, fcoeffs, tstamp, dir_out = tnet.get_network_parameters(net_name=net_name,experiment_name=net_name+'_PoA_demand')
        tNet = tnet.tNet(netFile=netFile, gFile=gFile, fcoeffs=fcoeffs)
    return tNet


v = []
x = []
i=0
for d in range(2,10):
    # System-centric
    tNet = read_net(net_name)
    g_per = tnet.perturbDemandConstant(tNet.g, d)
    tNet.set_g(g_per)
    tNet.build_supergraph(walk_multiplier=.001, identical_G=True)

    runtime = cars.solve_bush_CARSn(tNet, tNet.fcoeffs, n=5, rebalancing=False, bush=True, linear=False)
    #runtime = tNet.solveMSAsocial_supergraph(exogenous_G=False)
    cars.hist_flows(G=tNet.G_supergraph, G_exo=False)
    NLP_obj = tnet.get_totalTravelTime(tNet.G_supergraph, fcoeffs=tNet.fcoeffs, G_exogenous=False)
    del tNet
    # User-centric
    tNet = read_net(net_name)
    g_per = tnet.perturbDemandConstant(tNet.g,  d)
    tNet.set_g(g_per)
    tNet.build_supergraph(walk_multiplier=.001, identical_G=True)
    runtime = cars.solve_bush_CARSn(tNet, tNet.fcoeffs, n=5, rebalancing=False, bush=True, userCentric=True, linear=False)
    #runtime = tNet.solveMSA(exogenous_G=False)
    TAP_obj = tnet.get_totalTravelTime(tNet.G_supergraph, fcoeffs=tNet.fcoeffs, G_exogenous=False)

    del tNet

    poa = TAP_obj/NLP_obj
    v.append(poa)
    x.append(d*.2)
    print('d: ' +str(d*.2) + ', poa: '+ str(poa))

print(d)
print(v)

plt.plot(x,v, label='PoA')
plt.xlabel('demand')
plt.savefig('results/'+net_name  + '/_PoA_scan.pdf')
plt.show()