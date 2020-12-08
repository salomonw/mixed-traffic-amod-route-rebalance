import src.tnet as tnet
import src.CARS as cars
import numpy as np
import copy
from src.utils import *
import matplotlib as mpl
import matplotlib as mpl
from matplotlib import rc
import matplotlib.pyplot as plt
import experiments.build_NYC_subway_net as nyc
from datetime import datetime

#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'], 'size': 13})
rc('text', usetex=True)
plt.rc('axes', labelsize=13)
plt.rc('legend', fontsize=12)
plt.rc('xtick', labelsize=15)
plt.rc('ytick', labelsize=15)


def read_net(net_name):
    if net_name == 'NYC':
        tNet, tstamp, fcoeffs = nyc.build_NYC_net('data/net/NYC_M/', only_road=True)
    else:
        netFile, gFile, fcoeffs, tstamp, dir_out = tnet.get_network_parameters(net_name=net_name,
                                                                               experiment_name=net_name + '_n_variation')
        tNet = tnet.tNet(netFile=netFile, gFile=gFile, fcoeffs=fcoeffs)
    return tNet, fcoeffs

def solve_stackelberg_game(par):#netname, dir_out, penetration_rate, rebalancing=True, linear=False, n=5, theta_n=3, iterations=5):
    netname, dir_out, penetration_rate, rebalancing, linear, n, theta_n, iterations = par
    tNet_cavs, fcoeffs = read_net(netname)
    tNet_non_cavs = copy.deepcopy(tNet_cavs)

    g_cavs = tnet.perturbDemandConstant(tNet_cavs.g, constant=(penetration_rate))
    g_non_cavs = tnet.perturbDemandConstant(tNet_cavs.g, constant=(1 - penetration_rate))

    tNet_cavs.set_g(g_cavs)
    tNet_non_cavs.set_g(g_non_cavs)

    tNet_cavs.build_supergraph(walk_multiplier=1)#0.00001)
    tNet_non_cavs.build_supergraph(walk_multiplier=0.001)

    it = []
    for i in range(iterations):
        if i <=0:
            cars.solve_bush_CARSn(tnet=tNet_cavs, fcoeffs=fcoeffs, n=n, exogenous_G=False,
                                  rebalancing=True, bush=True, linear=linear,
                                  theta_n=theta_n)
        else:
            cars.solve_bush_CARSn(tnet=tNet_cavs, fcoeffs=fcoeffs, n=n, exogenous_G=tNet_non_cavs.G,
                                  rebalancing=True, bush=True, linear=linear,
                                  theta_n=theta_n)
        if i <=-10:
            cars.solve_bush_CARSn(tnet=tNet_non_cavs, fcoeffs=fcoeffs, n=n,
                                  rebalancing=False, bush=True, userCentric=True,
                                  exogenous_G=False, linear=linear,
                                  theta_n=theta_n)
        else:
            cars.solve_bush_CARSn(tnet=tNet_non_cavs, fcoeffs=fcoeffs, n=n,
                                  rebalancing=False, bush=True, userCentric=True,
                                  exogenous_G=tNet_cavs.G_supergraph, linear=linear,
                                  theta_n=theta_n)
        cars.supergraph2G(tNet_non_cavs)
        cars.G2supergraph(tNet_non_cavs)







        #cars.hist_flows(G=tNet_cavs.G_supergraph, G_exo=tNet_non_cavs.G_supergraph)
        '''We need to solve the issue of G vs superG'''
        cavscost_ = cars.get_totalTravelTime_without_Rebalancing(tNet_cavs, G_exogenous=tNet_non_cavs.G_supergraph)
        noncavscost_ = cars.get_totalTravelTime_without_Rebalancing(tNet_non_cavs, G_exogenous=tNet_cavs.G_supergraph)
        it.append((cavscost_+noncavscost_) / (
                tNet_cavs.totalDemand + tNet_non_cavs.totalDemand)*60)

        # PLOT ITERATION
        fig0 = plt.figure(figsize=(4, 3))
        ax0 = fig0.add_subplot()
        ax0.plot(it, linewidth=2, marker="o", label = '$\\gamma = $' + str(round(penetration_rate,1)), color='black' )
        plt.xlabel('Iteration number')
        plt.ylabel('$J$')
        ax0.grid(True)
        plt.tight_layout()
        ax0.legend()
        plt.xlim((0, 19))
        print(i)

    cavscost_= cars.get_totalTravelTime_without_Rebalancing(tNet_cavs, G_exogenous=tNet_non_cavs.G_supergraph)
    noncavscost_ = tnet.get_totalTravelTime(tNet_non_cavs.G, fcoeffs=fcoeffs, G_exogenous=tNet_cavs.G_supergraph)
    totalcost_ = (cavscost_+noncavscost_) / (tNet_cavs.totalDemand + tNet_non_cavs.totalDemand)*60
    cavscost_= cavscost_/tNet_cavs.totalDemand*60
    noncavscost_ = noncavscost_/ tNet_non_cavs.totalDemand*60


    cavsFlow = cars.get_amod_flow(tNet_cavs)
    nonCavsFlow = tnet.get_total_G_flow(tNet_non_cavs.G_supergraph)
    pedestrianFlow = cars.get_pedestrian_flow(tNet_cavs)
    rebalancingFlow = cars.get_rebalancing_flow(tNet_cavs)
    plt.tight_layout()
    mkdir_n('results/' + dir_out)
    if rebalancing == True:
        fig0.savefig('results/' + dir_out + '/iteration-rebalance-pR-+' + str(round(penetration_rate, 1)) + '.pdf')
    else:
        fig0.savefig('results/' + dir_out + '/iteration-pR-+' + str(round(penetration_rate, 1)) + '.pdf')
    print(penetration_rate)
    return cavscost_, noncavscost_, totalcost_, cavsFlow, nonCavsFlow, pedestrianFlow, rebalancingFlow

import multiprocessing as mp

def penRate(netname, dir_out, rebalancing=True, linear=False, n=5, theta_n=3):
    pool = mp.Pool(mp.cpu_count())
    cavsCost = []
    noncavsCost = []
    totCost = []
    cavsFlow = []
    nonCavsFlow = []
    pedestrianFlow = []
    rebalancingFlow = []
    fig0 = plt.figure(num=2, figsize=(4, 3))
    ax0 = fig0.add_subplot()
    '''
    for penetration_rate in np.linspace(0.01, 0.99, 11):
        tNet_cavs, fcoeffs = read_net(netname)
        tNet_non_cavs = copy.deepcopy(tNet_cavs)
        g_cavs = tnet.perturbDemandConstant(tNet_cavs.g, constant=(penetration_rate))
        g_non_cavs = tnet.perturbDemandConstant(tNet_cavs.g, constant=(1 - penetration_rate))
        tNet_cavs.set_g(g_cavs)
        tNet_non_cavs.set_g(g_non_cavs)
        tNet_cavs.build_supergraph(walk_multiplier=0.00001)#3.226)
        tNet_non_cavs.build_supergraph(walk_multiplier=0.00001)

        it = []
        for i in range(2):

            if i <=0:
                cars.solve_bush_CARSn(tnet=tNet_non_cavs, fcoeffs=fcoeffs, n=n,
                                      rebalancing=False, bush=True, userCentric=True,
                                      exogenous_G=False, linear=linear,
                                      theta_n=theta_n)
            else:
                cars.solve_bush_CARSn(tnet=tNet_non_cavs, fcoeffs=fcoeffs, n=n,
                                      rebalancing=False, bush=True, userCentric=True,
                                      exogenous_G=tNet_cavs.G_supergraph, linear=linear,
                                      theta_n=theta_n)

            cars.supergraph2G(tNet_non_cavs)
            cars.G2supergraph(tNet_non_cavs)

            if i <=-10:
                cars.solve_bush_CARSn(tnet=tNet_cavs, fcoeffs=fcoeffs, n=n, exogenous_G=False,
                                      rebalancing=True, bush=True, linear=linear,
                                      theta_n=theta_n)
            else:
                cars.solve_bush_CARSn(tnet=tNet_cavs, fcoeffs=fcoeffs, n=n, exogenous_G=tNet_non_cavs.G,
                                      rebalancing=True, bush=True, linear=linear,
                                      theta_n=theta_n)

            #cars.supergraph2G(tNet_non_cavs)


            #cars.supergraph2G(tNet_cavs)


            #cars.supergraph2G(tNet_cavs)
            #tNet_non_cavs.solveMSA(exogenous_G=tNet_cavs.G_supergraph)
            #cars.G2supergraph(tNet_non_cavs)

           #We need to solve the issue of G vs superG
            #cars.supergraph2G(tNet_non_cavs)
            #cavscost_ = tnet.get_totalTravelTime(tNet_cavs.G_supergraph, tNet_cavs.fcoeffs, tNet_non_cavs.G_supergraph)
            #noncavscost_= tnet.get_totalTravelTime(tNet_non_cavs.G_supergraph, tNet_non_cavs.fcoeffs, tNet_cavs.G_supergraph)
            cavscost_ = cars.get_totalTravelTime_without_Rebalancing(tNet_cavs, G_exogenous=tNet_non_cavs.G_supergraph)
            noncavscost_ = cars.get_totalTravelTime_without_Rebalancing(tNet_non_cavs, G_exogenous=tNet_cavs.G_supergraph)

            #noncavscost_ = tnet.get_totalTravelTime(tNet_non_cavs.G, fcoeffs=fcoeffs, G_exogenous=tNet_cavs.G_supergraph)
            #print(cavscost_/tNet_cavs.totalDemand*60)
            #print(noncavscost_/tNet_non_cavs.totalDemand*60)
            it.append((cavscost_+noncavscost_) / (
                    tNet_cavs.totalDemand + tNet_non_cavs.totalDemand)*60)

            if penetration_rate > -0.52 and penetration_rate < 5.53:
                fig0 = plt.figure(figsize=(4, 3))
                ax0 = fig0.add_subplot()
                ax0.plot(it, linewidth=2, marker="o", label = '$\\gamma = $' + str(round(penetration_rate,1)), color='black' )
                plt.xlabel('Iteration number')
                plt.ylabel('$J$')
                ax0.grid(True)
                plt.tight_layout()
                ax0.legend()
                plt.xlim((0, 19))

        #cars.supergraph2G(tNet_non_cavs)
        #cars.G2supergraph(tNet_non_cavs)

        #cavscost_ = tnet.get_totalTravelTime(tNet_cavs.G_supergraph, tNet_cavs.fcoeffs, tNet_non_cavs.G_supergraph)
        #noncavscost_= tnet.get_totalTravelTime(tNet_non_cavs.G_supergraph, tNet_non_cavs.fcoeffs, tNet_cavs.G_supergraph)
        cavscost_= cars.get_totalTravelTime_without_Rebalancing(tNet_cavs, G_exogenous=tNet_non_cavs.G_supergraph)
        #noncavscost_ = cars.get_totalTravelTime_without_Rebalancing(tNet_non_cavs, G_exogenous=tNet_cavs.G_supergraph)
        noncavscost_ = tnet.get_totalTravelTime(tNet_non_cavs.G, fcoeffs=fcoeffs, G_exogenous=tNet_cavs.G_supergraph)

        cavsCost.append(cavscost_/tNet_cavs.totalDemand*60)
        noncavsCost.append(noncavscost_/ tNet_non_cavs.totalDemand*60)
        totCost.append((cavscost_+noncavscost_) / (
                    tNet_cavs.totalDemand + tNet_non_cavs.totalDemand)*60)

        cavsFlow.append(cars.get_amod_flow(tNet_cavs))
        nonCavsFlow.append(tnet.get_total_G_flow(tNet_non_cavs.G_supergraph))
        pedestrianFlow.append(cars.get_pedestrian_flow(tNet_cavs))
        rebalancingFlow.append(cars.get_rebalancing_flow(tNet_cavs))

        mkdir_n('results/' + dir_out)
        if rebalancing == True:
            fig0.savefig('results/' + dir_out + '/iteration-rebalance-pR-+' + str(round(penetration_rate, 1)) + '.pdf')
        else:
            fig0.savefig('results/' + dir_out + '/iteration-rebalance-pR-+' + str(round(penetration_rate, 1)) + '.pdf')

        print(penetration_rate)
        print([cavsCost[-1], noncavsCost[-1], totCost[-1]])
        del tNet_cavs, tNet_non_cavs
    '''

    pars = [(netname, dir_out, penetration_rate, rebalancing, linear, n, theta_n, 5) for penetration_rate in np.linspace(0.001, 0.999, 11)]
    results = pool.map(solve_stackelberg_game, pars)
    i=0
    for penetration_rate in np.linspace(0.01, 0.99, 11):
        cavsCost.append(results[i][0])
        noncavsCost.append(results[i][1])
        totCost.append(results[i][2])
        cavsFlow.append(results[i][3])
        nonCavsFlow.append(results[i][4])
        pedestrianFlow.append(results[i][5])
        rebalancingFlow.append(results[i][6])
        i+=1

    return cavsCost, noncavsCost, totCost, cavsFlow, nonCavsFlow, pedestrianFlow, rebalancingFlow




# Inputa

#netname = 'NYC_Uber_small_1'
#netname = 'Braess1'
netname = 'EMA_mid'
#netname = 'EMA'
#netname = ' NYC'
#netname = 'Anaheim'
#netname = 'Sioux Falls'

rebalancing = True
n = 3


j = list(np.linspace(0.01, 0.99, 11))
alg = 'CARS'+str(n)
fig1 = plt.figure(num=2, figsize=(4, 2.8))
ax1 = fig1.add_subplot()
lstyle = ['-', '--', ':']
i=0


tstamp = datetime.now().strftime("%Y-%m-%d_%H:%M:%S")
dir_out = tstamp + '_' + 'penRate_'+netname

result = penRate(netname,dir_out=dir_out,rebalancing=rebalancing, n=n)
cavsCost, noncavsCost, totCost, cavsFlow, nonCavsFlow, pedestrianFlow, rebalancingFlow = result

ax1.plot(j, totCost, label=alg+' (Total)', color='green', linestyle=lstyle[i], linewidth=2, marker='o')
ax1.plot(j, cavsCost, label=alg+' (AMoDs)', color='blue', linestyle=lstyle[i], linewidth=2, marker="^")
ax1.plot(j, noncavsCost, label=alg+' (Private)', color='red', linestyle=lstyle[i], linewidth=2, marker='x')
ax1.legend()
plt.xlabel('Penetration Rate')
plt.ylabel('Avg. Travel Time (min)')
ax1.grid(True)
plt.tight_layout()
plt.xlim((0,1))
#plt.ylim((18,100))
plt.tight_layout()
mkdir_n('results/' + dir_out)
if rebalancing:
    fig1.savefig('results/' + dir_out +'/costs-rebalance.pdf')
else:
    fig1.savefig('results/' + dir_out +'/costs-no-rebalance.pdf')

fig2 = plt.figure(num=1, figsize=(6, 2.3))
ax2 = fig2.add_subplot()

fmt = mpl.ticker.StrMethodFormatter("{x}")
ax2.xaxis.set_major_formatter(fmt)
ax2.yaxis.set_major_formatter(fmt)

width = 0.9
x_name = [round(.1*i,1) for i in range(11)]
x = list(range(len(x_name)))
p1 = ax2.bar(x, nonCavsFlow, width)
p2 = ax2.bar(x, cavsFlow, width,
             bottom=nonCavsFlow)
p3 = ax2.bar(x, rebalancingFlow, width,
             bottom=[cavsFlow[i] + nonCavsFlow[i] for i in range(len(cavsFlow))])
p4 = ax2.bar(x, pedestrianFlow, width, bottom=[cavsFlow[i] + nonCavsFlow[i] + rebalancingFlow[i] for i in range(len(cavsFlow))])#, color='purple')
plt.ylabel('Flow')
plt.xlabel('Penetration rate', fontsize=10)
plt.xticks(x, x_name)#plt.legend((p1[0], p2[0], p3[0], p4[0]), ('Private Vehicle', 'AMoD', 'Rebalancing', 'Micromobility'), fontsize=12)
plt.legend((p1[0], p2[0], p3[0], p4[0]), ('Private Vehicle', 'AMoD', 'Rebalancing', 'Pedestrian'), fontsize=12)
ax2.tick_params(axis='both', which='major', labelsize=11)


plt.tight_layout()
mkdir_n('results/' + dir_out)

if rebalancing ==True:
    plt.savefig('results/' + dir_out +'/flows-rebalance.pdf')
else:
    plt.savefig('results/' + dir_out +'/flows-no-rebalance.pdf')
