import src.tnet as tnet
import src.CARS as cars
import numpy as np
import copy
from src.utils import *
import matplotlib as mpl
import matplotlib as mpl
from matplotlib import rc
import matplotlib.pyplot as plt
#rc('font',**{'family':'sans-serif','sans-serif':['He



def run_opt(penetration_rate, netFile, gFile, fcoeffs, alg='CARS3', rebalancing=True, xa=1):
    tNet_cavs = tnet.tNet(netFile=netFile, gFile=gFile, fcoeffs=fcoeffs)
    tNet_non_cavs = copy.deepcopy(tNet_cavs)
    g_cavs = tnet.perturbDemandConstant(tNet_cavs.g, constant=(penetration_rate)*2.5)
    g_non_cavs = tnet.perturbDemandConstant(tNet_cavs.g, constant=(1 - penetration_rate)*2.5)
    tNet_cavs.set_g(g_cavs)
    tNet_non_cavs.set_g(g_non_cavs)
    tNet_cavs.build_supergraph(walk_multiplier=0.00001)#3.226)
    tNet_non_cavs.build_supergraph(walk_multiplier=0.00001)
    fig0 = plt.figure(num=2, figsize=(4, 3))
    fig0.add_subplot()
    it = []
    for i in range(6):
        j=i
        if alg == 'CARS':
            if j == 0 :
                cars.solve_CARS(tNet_cavs, fcoeffs=fcoeffs, exogenous_G=False , xa=xa, rebalancing=rebalancing)
            else:
                cars.solve_CARS(tNet_cavs, fcoeffs=fcoeffs, exogenous_G=tNet_non_cavs.G_supergraph,  xa=xa,
                                rebalancing=rebalancing)
            cars.supergraph2G(tNet_cavs)
        elif alg == 'CARS3':
            if j==0:
                a, b, od_flows =cars.solve_CARS2(tNet_cavs, fcoeffs=fcoeffs, exogenous_G=False,  rebalancing=rebalancing)
            else:
                a, b, od_flows = cars.solve_CARS2(tNet_cavs, fcoeffs=fcoeffs, exogenous_G=tNet_non_cavs.G_supergraph,
                                 rebalancing=rebalancing)
        elif alg == 'disjoint':
            if j == 0:
                #cars.solve_social_Julia(tNet, exogenous_G=0)
                cars.solve_social_Julia(tNet_cavs, exogenous_G=False)
                if rebalancing != False:
                    cars.supergraph2G(tNet_cavs)
                    cars.solve_rebalancing(tNet_cavs, exogenous_G=tNet_non_cavs.G_supergraph)
            else:
                cars.solve_social_Julia(tNet_cavs, exogenous_G=tNet_non_cavs.G_supergraph)
                if rebalancing != False:
                    cars.supergraph2G(tNet_cavs)
                    cars.solve_rebalancing(tNet_cavs, exogenous_G=tNet_non_cavs.G_supergraph)
        else:
            print('please select an algorithm (CARS, CARS2 or disjoint) !!!!!!!!!')
            break

        tNet_non_cavs.solveMSA(exogenous_G=tNet_cavs.G_supergraph)
        cars.G2supergraph(tNet_non_cavs)

        totalCost = cars.eval_obj_funct(tNet_cavs, G_exogenous=tNet_non_cavs.G)
        it.append(totalCost)

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

    cavscost_= cars.get_totalTravelTime_without_Rebalancing(tNet_cavs, G_exogenous=tNet_non_cavs.G_supergraph)
    noncavscost_ = tnet.get_totalTravelTime(tNet_non_cavs.G_supergraph, fcoeffs=fcoeffs, G_exogenous=tNet_cavs.G_supergraph)

    mkdir_n('results/' + dir_out)
    if rebalancing == True:
        fig0.savefig('results/' + dir_out + '/iteration-rebalance-pR-+' + str(round(penetration_rate, 1)) + '.pdf')
    else:
        fig0.savefig('results/' + dir_out + '/iteration-rebalance-pR-+' + str(round(penetration_rate, 1)) + '.pdf')
    return cavscost_, noncavscost_, tNet_cavs, od_flows



import src.routeFinder as routeFinder

netFile, gFile, fcoeffs, tstamp, dir_out = tnet.get_network_parameters('EMA', experiment_name='EMA_penRate_comparison-'+'REB')
cavscost_, noncavscost_, tNet_cavs, od_flows = run_opt(0.5, netFile, gFile, fcoeffs, alg='CARS3', rebalancing=True, xa=1)
tNet_cavs.set_supergraph_tk(fcoeffs=fcoeffs, G_exogenous=False)
routes = routeFinder.routeFinder(tNet_cavs.G_supergraph, tNet_cavs.g, od_flows, eps=20)
print(routes)

#print(od_flows)
plt.show()

