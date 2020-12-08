import src.tnet as tnet
import src.CARS as cars
import numpy as np
import copy
from src.utils import *
import matplotlib as mpl
from matplotlib import rc
import matplotlib.pyplot as plt
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

#netFile, gFile, fcoeffs, tstamp, dir_out = tnet.get_network_parameters('EMA', experiment_name='EMA_penRate_comparison-'+'REB')

netFile, gFile, fcoeffs, tstamp, dir_out = tnet.get_network_parameters('NYC_Uber_small_1', experiment_name='NYC_Uber_small_1_penRate_comparison-REB')
#netFile, gFile, fcoeffs, tstamp, dir_out = tnet.get_network_parameters('Anaheim', experiment_name='Anaheim_test_CARSn')
#netFile, gFile, fcoeffs, tstamp, dir_out = tnet.get_network_parameters('Anaheim', experiment_name='Anaheim_test_CARSn')


def runtime_algs():
    probs = {'LP':[0,1,2,3,4,5], 'QP':[1,2]}

    runtimes = {}
    for prob in probs.keys():
        runtimes[prob] = {}
        for alg in probs[prob]:
            runtimes[prob][alg] = {}
            runtimes[prob][alg]['runtime'] = []
            runtimes[prob][alg]['obj'] = []

    for i in range(6):
        theta, a = cars.get_approx_fun(fcoeffs=fcoeffs, nlines=i+2, range_=[0, 3], plot=False)
        theta.append(3)
        print('\n')
        for prob, algs in probs.items():
            for alg in algs:
                tNet = tnet.tNet(netFile=netFile, gFile=gFile, fcoeffs=fcoeffs)
                tNet.build_supergraph(walk_multiplier=1)
                q = [(u, v) for (u, v, d) in tNet.G_supergraph.edges(data=True) if d['type'] == 'p']
                q = [(u, v) for (u, v, d) in tNet.G_supergraph.edges(data=True) if d['type'] == 'f']
                g_per = tnet.perturbDemandConstant(tNet.g, 1)
                tNet.set_g(g_per)

                if prob == 'LP':
                    linear = True
                else:
                    linear = False
                tNet, runtime, od_flows = cars.solve_CARSn(tNet, fcoeffs=fcoeffs,
                                                           n=i+2, exogenous_G=False,
                                                           rebalancing=True, linear=linear,
                                                           theta=theta, a=a, method=alg)

                obj = tnet.get_totalTravelTime(tNet.G_supergraph, fcoeffs)

                runtimes[prob][alg]['runtime'].append(runtime)
                runtimes[prob][alg]['obj'].append(obj)

                print('Prob:' +prob+ '\t | \t Alg: ' + str(alg) + '\t | \t n: ' + str(i + 2) + '\t | \t runtime: ' + str(round(runtime,3))+ '\t | \t obj: ' + str(round(obj,0)))

            del tNet
    ns = []
    for i in range(6):
        ns.append(i+2)
    return runtimes, ns




runtimes, ns = runtime_algs()
algs_dict = {0:'primal spmx', 1:'dual spmx', 2:'barrier', 3:'concurrent', \
             4:'det concurrent', 5:'det concurrent spmx'}

fig, ax = plt.subplots(4, figsize=(6,8))
for prob,v in runtimes.items():
    if prob == 'LP':
        for alg,j in v.items():
            ax[0].plot(ns, j['runtime'], label=algs_dict[alg])
            if alg==1:
                ax[2].plot(ns, j['runtime'], label=prob+'-'+algs_dict[alg])
        ax[3].plot(ns, j['obj'], label=prob)
    else:
        for alg,j in v.items():
            ax[1].plot(ns, j['runtime'], label=algs_dict[alg])
            if alg==2:
                ax[2].plot(ns, j['runtime'], label=prob+'-'+algs_dict[alg])
        ax[3].plot(ns, j['obj'], label=prob)

for i in ax:
    i.grid(True)
    i.legend()
    i.set_xlim([ns[0], ns[-1]])
ax[0].set_ylabel('LP runtime')
ax[1].set_ylabel('QP runtime')
ax[2].set_ylabel('runtime')
ax[3].set_ylabel('Obj Function')
plt.tight_layout()

plt.show()

