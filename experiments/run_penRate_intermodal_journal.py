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
import  pandas as pd
import multiprocessing as mp
import sys

#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'], 'size': 13})
rc('text', usetex=False)
plt.rc('axes', labelsize=13)
plt.rc('legend', fontsize=12)
plt.rc('xtick', labelsize=15)
plt.rc('ytick', labelsize=15)


def read_net(net_name):
    if net_name == 'NYC':
        tNet, tstamp, fcoeffs = nyc.build_NYC_net('data/net/NYC/', only_road=True)
    else:
        netFile, gFile, fcoeffs, tstamp, dir_out = tnet.get_network_parameters(net_name=net_name,
                                                                               experiment_name=net_name + '_n_variation')
        tNet = tnet.tNet(netFile=netFile, gFile=gFile, fcoeffs=fcoeffs)
    return tNet, fcoeffs

def get_pwfunction(fcoeffs, n, theta_n,  userCentric=False):
    fc = fcoeffs.copy()
    if userCentric:
        fc = cars.UC_fcoeffs(fc)
    theta, a, rms  = cars.get_approx_fun(fcoeffs=fc, nlines=n, range_=[0,theta_n], plot=False)
    return theta, a

def solve_stackelberg_game(par):
    netname, dir_out, penetration_rate, rebalancing, linear, n, theta_n, demand_multiplier,  iterations = par
    tNet_cavs, fcoeffs = read_net(netname)
    tNet_non_cavs = copy.deepcopy(tNet_cavs)
    g_cavs = tnet.perturbDemandConstant(tNet_cavs.g, constant=(penetration_rate*demand_multiplier))
    g_non_cavs = tnet.perturbDemandConstant(tNet_cavs.g, constant=(1 - penetration_rate)*demand_multiplier)

    tNet_cavs.set_g(g_cavs)
    tNet_non_cavs.set_g(g_non_cavs)

    tNet_cavs.build_supergraph(identical_G=True)
    tNet_non_cavs.build_supergraph(identical_G=True)

    theta_cavs, a_cavs = get_pwfunction(fcoeffs, n, theta_n, userCentric=False)
    theta_non_cavs, a_non_cavs = get_pwfunction(fcoeffs, n, theta_n, userCentric=True)

    #'''
    # add biking network
    if "B" in str(sys.argv[3]):
        tNet_cavs.build_layer(one_way=True, avg_speed=6,capacity=99999, symb="b", identical_G=False)
        tNet_non_cavs.build_layer(one_way=True, avg_speed=0.0001, capacity=99999, symb="b", identical_G=False)
    # add pedestrian network
    if "P" in str(sys.argv[3]):
        tNet_cavs.build_layer(one_way=False, avg_speed=3.1, capacity=99999, symb="'", identical_G=False)
        tNet_non_cavs.build_layer(one_way=False, avg_speed=0.0001, capacity=9999, symb="'", identical_G=False)
    # add subway for nyc
    if netname == 'NYC' and "S" in str(sys.argv[3]):
        layer = tnet.readNetFile(netFile='data/net/NYC/NYC_M_Subway_net.txt')
        tNet_cavs.add_layer(layer=layer, layer_symb='s')
        tNet_non_cavs.add_layer(layer=layer, layer_symb='s', speed=0.001)
  
    it = []
    for i in range(iterations):
       # '''
        if i <= 0:
            cars.solve_bush_CARSn(tnet=tNet_non_cavs, fcoeffs=fcoeffs, n=n, exogenous_G=False,
	                          rebalancing=rebalancing, bush=True, linear=linear,
				  theta_n=theta_n, theta=theta_non_cavs, a=a_non_cavs, od_flows_flag=False, 
				  userCentric=True)
        else:
            cars.solve_bush_CARSn(tnet=tNet_non_cavs, fcoeffs=fcoeffs, n=n, exogenous_G=tNet_cavs.G_supergraph, 
				  rebalancing=rebalancing, bush=True, linear=linear, 
				  theta_n=theta_n, theta=theta_non_cavs, a=a_non_cavs, od_flows_flag=False, 
				  userCentric=True)
       # '''
        if i <=-10:
            cars.solve_bush_CARSn(tnet=tNet_cavs, fcoeffs=fcoeffs, n=n, exogenous_G=False,
                                  rebalancing=rebalancing, bush=True, linear=linear,
                                  theta_n=theta_n, theta=theta_cavs, a=a_cavs, od_flows_flag=False)
        else:
            cars.solve_bush_CARSn(tnet=tNet_cavs, fcoeffs=fcoeffs, n=n, exogenous_G=tNet_non_cavs.G_supergraph,
                                  rebalancing=rebalancing, bush=True, linear=linear,
                                  theta_n=theta_n, theta=theta_cavs, a=a_cavs, od_flows_flag=False)
        ''' 
        if i <=-10:
            #tNet_non_cavs.solveMSA(exogenous_G=False)
            cars.solve_bush_CARSn(tnet=tNet_non_cavs, fcoeffs=fcoeffs, n=n,
                                  rebalancing=False, bush=True, userCentric=True,
                                  exogenous_G=False, linear=linear,
                                  theta_n=theta_n, theta=theta_non_cavs, a=a_non_cavs, od_flows_flag=False)
        else:
            #tNet_non_cavs.solveMSA(exogenous_G=tNet_cavs.G_supergraph)
            cars.solve_bush_CARSn(tnet=tNet_non_cavs, fcoeffs=fcoeffs, n=n,
                                  rebalancing=False, bush=True, userCentric=True,
                                  exogenous_G=tNet_cavs.G_supergraph, linear=linear,
                                  theta_n=theta_n, theta=theta_non_cavs, a=a_non_cavs, od_flows_flag=False)
        '''  
        #cars.G2supergraph(tNet_non_cavs)
        #cars.supergraph2G(tNet_non_cavs)

        cavscost_ = cars.get_totalTravelTime_without_Rebalancing(tNet_cavs, G_exogenous=tNet_non_cavs.G_supergraph)
        noncavscost_ = cars.get_totalTravelTime_without_Rebalancing(tNet_non_cavs, G_exogenous=tNet_cavs.G_supergraph)
        it.append((cavscost_+noncavscost_) / (tNet_cavs.totalDemand + tNet_non_cavs.totalDemand)*60)


        # PLOT ITERATION
        fig0 = plt.figure(figsize=(4, 3))
        ax0 = fig0.add_subplot()
        ax0.plot(it, linewidth=2, marker="o", label = '$\\gamma = $' + str(round(penetration_rate,1)), color='black' )
        plt.xlabel('Iteration number')
        plt.ylabel('$J$')
        ax0.grid(True)
        plt.tight_layout()
        ax0.legend()
        plt.xlim((0, iterations-1))
        print(i)

        cars.hist_flows(G=tNet_cavs.G_supergraph, G_exo=tNet_non_cavs.G_supergraph)

    cavscost_= cars.get_totalTravelTime_without_Rebalancing(tNet_cavs, G_exogenous=tNet_non_cavs.G_supergraph)
    noncavscost_ = cars.get_totalTravelTime_without_Rebalancing(tNet_non_cavs, G_exogenous=tNet_cavs.G_supergraph)
    totalcost_ = (cavscost_+noncavscost_) / (tNet_cavs.totalDemand + tNet_non_cavs.totalDemand)*60
    cavscost_= cavscost_/tNet_cavs.totalDemand*60
    noncavscost_ = noncavscost_/ tNet_non_cavs.totalDemand*60



    cavsFlow = cars.get_amod_flow(tNet_cavs)
    rebalancingFlow = cars.get_rebalancing_flow(tNet_cavs)
    nonCavsFlow = tnet.get_total_G_flow(tNet_non_cavs.G_supergraph)
    pedestrianFlow = cars.get_layer_flow(tNet_cavs, symb="'")
    bikeFlow = cars.get_layer_flow(tNet_cavs, symb="b")
    subwayFlow = 0
    if netname == 'NYC':
        subwayFlow = cars.get_layer_flow(tNet_cavs, symb="s")
    print(subwayFlow)
    identical_G=False
    plt.tight_layout()
    mkdir_n('results/' + dir_out)
    if rebalancing == True:
        fig0.savefig('results/' + dir_out + '/iteration-rebalance-pR-+' + str(round(penetration_rate, 1)) + '.pdf')
    else:
        fig0.savefig('results/' + dir_out + '/iteration-pR-+' + str(round(penetration_rate, 1)) + '.pdf')
    print(penetration_rate)
    return cavscost_, noncavscost_, totalcost_, cavsFlow, nonCavsFlow, pedestrianFlow, rebalancingFlow, bikeFlow, subwayFlow

def penRate(netname, dir_out, rebalancing=True, linear=False, n=5, theta_n=3, demand_multiplier=1,n_iter = 5,  parallel=True):

    j  = [0.01, 0.1, 0.2, 0.3 , 0.4, 0.5 ,0.6, 0.7 , 0.8 , 0.9, 0.99]
    cavsCost = []
    noncavsCost = []
    totCost = []
    cavsFlow = []
    nonCavsFlow = []
    pedestrianFlow = []
    rebalancingFlow = []
    bikeFlow = []
    subwayFlow = []
    if parallel:
        pool = mp.Pool(12)
        pars = [(netname, dir_out, penetration_rate, rebalancing, linear, n, theta_n, demand_multiplier, n_iter) for penetration_rate in j]
        results = pool.map(solve_stackelberg_game, pars)
    else:
        results = []
        for penetration_rate in j:
            par = (netname, dir_out, penetration_rate, rebalancing, linear, n, theta_n, demand_multiplier, n_iter)
            results.append(solve_stackelberg_game(par))

    i=0
    for penetration_rate in j:
        cavsCost.append(results[i][0])
        noncavsCost.append(results[i][1])
        totCost.append(results[i][2])
        cavsFlow.append(results[i][3])
        nonCavsFlow.append(results[i][4])
        pedestrianFlow.append(results[i][5])
        rebalancingFlow.append(results[i][6])
        bikeFlow.append(results[i][7])
        subwayFlow.append(results[i][8])
        i+=1

    return cavsCost, noncavsCost, totCost, cavsFlow, nonCavsFlow, pedestrianFlow, rebalancingFlow, bikeFlow, subwayFlow

def plot_penRate(result, dir_out):
    cavsCost, noncavsCost, totCost, cavsFlow, nonCavsFlow, pedestrianFlow, rebalancingFlow, bikeFlow, subwayFlow = result

    fig1 = plt.figure(num=2, figsize=(4, 2.8))
    ax1 = fig1.add_subplot()
    lstyle = ['-', '--', ':']
    i = 0
    ax1.plot(j, totCost, label=alg + ' (Total)', color='green', linestyle=lstyle[i], linewidth=2, marker='o')
    ax1.plot(j, cavsCost, label=alg + ' (AMoDs)', color='blue', linestyle=lstyle[i], linewidth=2, marker="^")
    ax1.plot(j, noncavsCost, label=alg + ' (Private)', color='red', linestyle=lstyle[i], linewidth=2, marker='x')
    ax1.legend()
    plt.xlabel('Penetration Rate')
    plt.ylabel('Avg. Travel Time (min)')
    ax1.grid(True)
    plt.tight_layout()
    plt.xlim((0, 1))
    # plt.ylim((18,100))
    plt.tight_layout()
    mkdir_n('results/' + dir_out)
    if rebalancing:
        fig1.savefig('results/' + dir_out + '/costs-rebalance.pdf')
    else:
        fig1.savefig('results/' + dir_out + '/costs-no-rebalance.pdf')

    fig2 = plt.figure(num=1, figsize=(6, 2.3))
    ax2 = fig2.add_subplot()

    fmt = mpl.ticker.StrMethodFormatter("{x}")
    ax2.xaxis.set_major_formatter(fmt)
    ax2.yaxis.set_major_formatter(fmt)

    width = 0.9
    x_name = [round(.1 * i, 1) for i in range(11)]
    x = list(range(len(x_name)))
    p1 = ax2.bar(x, nonCavsFlow, width)
    p2 = ax2.bar(x, cavsFlow, width,
                 bottom=nonCavsFlow)
    p3 = ax2.bar(x, rebalancingFlow, width,
                 bottom=[cavsFlow[i] + nonCavsFlow[i] for i in range(len(cavsFlow))])
    p4 = ax2.bar(x, pedestrianFlow, width,
                 bottom=[cavsFlow[i] + nonCavsFlow[i] + rebalancingFlow[i] for i in range(len(cavsFlow))])
    p5 = ax2.bar(x, bikeFlow, width,
                 bottom=[cavsFlow[i] + nonCavsFlow[i] + rebalancingFlow[i] + pedestrianFlow[i] for i in
                         range(len(cavsFlow))])
    p6 = ax2.bar(x, subwayFlow, width,
                 bottom=[cavsFlow[i] + nonCavsFlow[i] + rebalancingFlow[i] + pedestrianFlow[i] + bikeFlow[i] for i in
                         range(len(cavsFlow))])

    plt.ylabel('Flow')
    plt.xlabel('Penetration rate', fontsize=10)
    plt.xticks(x, x_name)  # plt.legend((p1[0], p2[0], p3[0], p4[0]), ('Private Vehicle', 'AMoD', 'Rebalancing', 'Micromobility'), fontsize=12)
    plt.legend((p1[0], p2[0], p3[0], p4[0], p5[0], p6[0]), ('Private Vehicle', 'AMoD', 'Rebalancing', 'Pedestrian', 'Bike', 'Subway'),
               fontsize=12)
    ax2.tick_params(axis='both', which='major', labelsize=11)

    plt.tight_layout()
    mkdir_n('results/' + dir_out)

    if rebalancing == True:
        plt.savefig('results/' + dir_out + '/flows-rebalance.pdf')
    else:
        plt.savefig('results/' + dir_out + '/flows-no-rebalance.pdf')

def save_results(results, dir_out):
    labels = ['cavsCost', 'noncavsCost', 'totCost', 'cavsFlow', 'nonCavsFlow', 'pedestrianFlow', 'rebalancingFlow', 'bikeFlow', 'subwayFlow']
    df = pd.DataFrame(list(map(list, zip(*results))), columns=labels)
    df.to_csv('results/' + dir_out + '/results.csv')



# Inputs
netname = str(sys.argv[1])
demand_multiplier = float(sys.argv[2])
modes = str(sys.argv[3])
rebalancing = True
n = 5
n_iter = 5
theta_n = 1.5
linear = False
parallel = True

if netname == "NYC":
    parallel = True
j = [0, 0.1, 0.2, 0.3 , 0.4, 0.5 ,0.6, 0.7 , 0.8 , 0.9, 1]
alg = 'CARS'+str(n)
tstamp = datetime.now().strftime("%Y-%m-%d_%H:%M:%S")
dir_out = tstamp + '_' + 'penRate_'+netname+"_"+str(demand_multiplier)+modes
result = penRate(netname,
		dir_out=dir_out,
		rebalancing=rebalancing, 
		n=n, 
		theta_n=theta_n, 
		linear=linear,
		demand_multiplier=demand_multiplier, 
		n_iter=n_iter, 
		parallel=parallel
		)
save_results(result, dir_out)
plot_penRate(result, dir_out)
#plt.show()
with open('results/' + dir_out + "/parameters.txt", "w") as text_file:
            print("n: "+str(n)+"\ntheta_n: "+str(theta_n)+"\nlinear: "+str(linear)+"\ng_multi: "+str(demand_multiplier), file=text_file)
