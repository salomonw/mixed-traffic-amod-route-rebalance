from src.utils import *
import src.tnet as tnet
import src.CARS as cars
import matplotlib as mpl
import numpy as np
import json
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib import rc
import matplotlib.pyplot as plt
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)


#netFile, gFile, fcoeffs, tstamp, dir_out = tnet.get_network_parameters('EMA', experiment_name='EMA_penRate_comparison')

#netFile, gFile, fcoeffs, tstamp, dir_out = tnet.get_network_parameters('NYC_Uber_small', experiment_name='NYC_Uber_small_penRate_comparison')
xa = 0.8

def computation_time(netFile, gFile, fcoeffs, g=[0.8,1.2], alg='CARS3', n=30, runtimesDir={}):
    comp_time = []
    obj_v = []
    gmutli = [np.random.uniform(g[0],g[1]) for i in range(n)]
    for i in range(n):
        tNet = tnet.tNet(netFile=netFile, gFile=gFile, fcoeffs=fcoeffs)
        tNet.build_supergraph(walk_multiplier=999)
        [(u, v) for (u, v, d) in tNet.G_supergraph.edges(data=True) if d['type'] == 'p']
        [(u, v) for (u, v, d) in tNet.G_supergraph.edges(data=True) if d['type'] == 'f']
        g_per = tnet.perturbDemandConstant(tNet.g, gmutli[i])
        tNet.set_g(g_per)

        if alg=='CARS':
            tNet, t = cars.solve_CARS(tNet, exogenous_G=0, fcoeffs=fcoeffs, xa=1.2, rebalancing=True)
            obj = cars.eval_obj_funct(tNet, G_exogenous=False)
            runtimes['CARS'].append(t)
        elif alg == 'CARS3':
            tNet, t = cars.solve_CARS2(tNet, exogenous_G=0, fcoeffs=fcoeffs, rebalancing=True)
            obj = cars.eval_obj_funct(tNet, G_exogenous=False)
            runtimes['CARS3'].append(t)
        elif alg=='disjoint':
            t_NLP = cars.solve_social_Julia(tNet, exogenous_G=0)
            cars.supergraph2G(tNet)
            t_Reb = cars.solve_rebalancing(tNet, exogenous_G=0)
            t = t_NLP + t_Reb
            obj = cars.eval_obj_funct(tNet, G_exogenous=False)
            runtimes['NLP'].append(t_NLP)
            runtimes['Reb'].append(t_Reb)
            runtimes['Disjoint'].append(t)
        else:
            print('invalid algorithm name')
        comp_time.append(t)
        obj_v.append(obj)
    return comp_time, obj_v, tNet.totalDemand


print('Runtimes Experiment\n------------------------------------------------------------')
print('Net \t Alg \t n \t Mean (min) \t Var \t Obj (J/sum(g))')
print('------------------------------------------------------------')

n = 30
runtimes = {}
for net in ['EMA', 'NYC']:
    # Set network
    if net =='EMA':
        netFile, gFile, fcoeffs, tstamp, dir_out = tnet.get_network_parameters('EMA',experiment_name='EMA_runtime')
    elif net == 'NYC':
        netFile, gFile, fcoeffs, tstamp, dir_out = tnet.get_network_parameters('NYC_Uber_small_1',
                                                                               experiment_name='NYC_runtime')
    elif net == "Anaheim":
        netFile, gFile, fcoeffs, tstamp, dir_out = tnet.get_network_parameters('Anaheim', experiment_name='Anaheim_runtime')
    else:
        print('net acronym wrong')
    # set random demand between two values


    # solve the problem for different algoritms
    for alg in ['CARS', 'CARS3', 'disjoint']:

        if alg == 'disjoint':
            runtimes['NLP'] = []
            runtimes['Reb'] = []
            runtimes['Disjoint'] = []
        else:
            runtimes[alg] = []

        t, obj, totalDemand = computation_time(netFile, gFile, fcoeffs, alg=alg, n=n, runtimesDir=runtimes)
        print(net + '\t' + alg + '\t' + str(n) +'\t' + str(round(np.mean(t),6)) +'\t' + str(round(np.var(t),8)) +'\t' + str(round(np.mean(obj),3)))
print('------------------------------------------------------------')


# save results
mkdir_n('results/' + tstamp + '_runtimes_experiment/')
js = json.dumps(runtimes)
f = open('results/' + tstamp + '_runtimes_experiment/runtimes.txt', "w")
f.write(js)
f.close()

print('Run times saved to: \t' + 'results/' + 'results/' + tstamp + '/runtimes_experiment')
