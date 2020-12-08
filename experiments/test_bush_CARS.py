import src.tnet as tnet
import src.CARS as cars
from matplotlib import rc
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import pandas as pd
import experiments.build_NYC_subway_net as nyc
import random as random
from src.utils import mkdir_n

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#rc('text', usetex=True)


def run_n_variation(net_name, g_mult):
    n = [2+i for i in range(7)]
    CARS = {}
    d = []
    TAP_obj = 1
    NLP_obj = 1

    if net_name=='NYC':
        tNet, tstamp, fcoeffs = nyc.build_NYC_net('data/net/NYC_M/', only_road=True)
    else:
        netFile, gFile, fcoeffs, tstamp, dir_out = tnet.get_network_parameters(net_name=net_name,
                                                                    experiment_name=net_name+'_n_variation')
    #'''
    if net_name != 'NYC':
        tNet = tnet.tNet(netFile=netFile, gFile=gFile, fcoeffs=fcoeffs)
    else:
        tNet, tstamp, fcoeffs = nyc.build_NYC_net('data/net/NYC_M/', only_road=True)
    g_per = tnet.perturbDemandConstant(tNet.g, g_mult)
    tNet.set_g(g_per)
    tNet.build_supergraph(walk_multiplier=.001, identical_G=True)
    runtime, RG = cars.solveMSAsocialCARS(tNet, exogenous_G=False)
    NLP_obj = tnet.get_totalTravelTime(tNet.G_supergraph, fcoeffs=fcoeffs, G_exogenous=False)
    d0 = {'Net': net_name, 'A': tNet.nLinks, 'V': tNet.nNodes, 'W': tNet.nOD, 'type': 'system-centric', 'obj':NLP_obj, 't':runtime, 'RG':RG}
    d.append(d0)

    del tNet
    if net_name != 'NYC':
        tNet = tnet.tNet(netFile=netFile, gFile=gFile, fcoeffs=fcoeffs)
    else:
        tNet, tstamp, fcoeffs = nyc.build_NYC_net('data/net/NYC_M/', only_road=True)
    g_per = tnet.perturbDemandConstant(tNet.g, g_mult)
    tNet.set_g(g_per)
    tNet.build_supergraph(walk_multiplier=.001, identical_G=True)
    runtime, RG = tNet.solveMSA(exogenous_G=False)
    TAP_obj = tnet.get_totalTravelTime(tNet.G, fcoeffs=fcoeffs, G_exogenous=False)
    d0 = {'Net': net_name, 'A': tNet.nLinks, 'V': tNet.nNodes, 'W': tNet.nOD, 'type': 'user-centric', 'obj':TAP_obj, 't':runtime, 'RG':RG}
    d.append(d0)
    #'''
    tap_obj =[]
    for i in n:
        CARS[i] = {}
        tap_obj.append((TAP_obj/NLP_obj-1)*100)
        for linear in [True, False]:

            if net_name != 'NYC':
                tNet = tnet.tNet(netFile=netFile, gFile=gFile, fcoeffs=fcoeffs)
            else:
                tNet, tstamp, fcoeffs = nyc.build_NYC_net('data/net/NYC_M/')
            g_per = tnet.perturbDemandConstant(tNet.g, g_mult)
            tNet.set_g(g_per)
            tNet.build_supergraph(walk_multiplier=.001, identical_G=True)
            tNet, runtime = cars.solve_bush_CARSn(tNet, fcoeffs=fcoeffs, n=i, exogenous_G=False, rebalancing=False, linear=linear,  bush=True, theta_n=5, userCentric=False, od_flows_flag=False)
            CARS2obj = tnet.get_totalTravelTime(tNet.G_supergraph, fcoeffs)
            CARS[i][linear] = (CARS2obj/NLP_obj-1)*100
            if linear == True:
                tp = 'LP'
            else:
                tp = 'QP'
            d0 = {'Net': net_name, 'A': tNet.nLinks, 'V': tNet.nNodes, 'W': tNet.nOD, 'type': 'CARS'+str(i)+'-'+tp, 'obj':CARS2obj, 't':runtime, 'n':i, 'g':g_mult}
            d.append(d0)
            print(pd.DataFrame(d))
            del tNet

    fig, ax = plt.subplots(figsize=(5,2))
    ax.plot(n, [v[True] for k,v in CARS.items()], label = 'LP', marker=".")
    ax.plot(n, [v[False] for k,v in CARS.items()], label = 'QP', marker=">")
    ax.plot(n, tap_obj, label = 'user-centric', marker="x", linestyle="--")
    ax.set_xlabel('n')
    #ax.set_ylabel('% deviation from NLP')
    ax.set_xlim([n[0], n[-1]])
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.legend(framealpha=1)
    ax.grid(True)
    plt.tight_layout()
    plt.savefig(outdir+"/"+net_name+'_'+str(g_mult)+'.pdf')
    #plt.show()
    return d


#nets = ['NYC_Uber_small_1']#, 'NYC_Uber_small_1','EMA_mid', 'NYC']
#nets = ['Sioux Falls', 'Anaheim', 'EMA_mid', 'NYC']
#nets = ['EMA_mid', 'NYC']
nets = ['NYC']
#nets = ['EMA']
#nets = ['Sioux Falls', 'Anaheim', 'ChicagoSketch']
g_mult = [4]
d_big = []
outdir = 'results/variation/'
mkdir_n(outdir)

for net in nets:
    for g in g_mult:
        d = run_n_variation(net, g)
        d_big.extend(d)
        df = pd.DataFrame(d)
        df.to_csv(outdir+'/n_variation_'+net+'_g_'+str(g)+'_user.csv')
df_big = pd.DataFrame(d_big)

