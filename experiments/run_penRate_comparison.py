import src.tnet as tnet
import src.CARS as cars
import numpy as np
import copy
from src.utils import *
import matplotlib as mpl
import matplotlib as mpl
from matplotlib import rc
import matplotlib.pyplot as plt
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'], 'size': 13})
rc('text', usetex=True)
plt.rc('axes', labelsize=13)
plt.rc('legend', fontsize=12)
plt.rc('xtick', labelsize=15)
plt.rc('ytick', labelsize=15)


def penRate(netFile, gFile, fcoeffs, alg='CARS3', rebalancing=True, xa=1):
    cavsCost = []
    noncavsCost = []
    totCost = []
    cavsFlow = []
    nonCavsFlow = []
    pedestrianFlow = []
    rebalancingFlow = []
    fig0 = plt.figure(num=2, figsize=(4, 3))
    ax0 = fig0.add_subplot()
    for penetration_rate in np.linspace(0.01, 0.99, 11):
        tNet_cavs = tnet.tNet(netFile=netFile, gFile=gFile, fcoeffs=fcoeffs)
        tNet_non_cavs = copy.deepcopy(tNet_cavs)
        g_cavs = tnet.perturbDemandConstant(tNet_cavs.g, constant=(penetration_rate)*2.5)
        g_non_cavs = tnet.perturbDemandConstant(tNet_cavs.g, constant=(1 - penetration_rate)*2.5)
        tNet_cavs.set_g(g_cavs)
        tNet_non_cavs.set_g(g_non_cavs)
        tNet_cavs.build_supergraph(walk_multiplier=0.00001)#3.226)
        tNet_non_cavs.build_supergraph(walk_multiplier=0.00001)

        it = []
        for i in range(6):

            j=i
            '''
            if i == 0 :
                tNet_non_cavs.solveMSA(exogenous_G=0)
                cars.G2supergraph(tNet_non_cavs)
                j+=1
            else:
                tNet_non_cavs.solveMSA(exogenous_G=tNet_cavs.G_supergraph)
                cars.G2supergraph(tNet_non_cavs)
             '''


            if alg == 'CARS':
                if j == 0 :
                    cars.solve_CARS(tNet_cavs, fcoeffs=fcoeffs, exogenous_G=False , xa=xa, rebalancing=rebalancing)
                else:
                    cars.solve_CARS(tNet_cavs, fcoeffs=fcoeffs, exogenous_G=tNet_non_cavs.G_supergraph,  xa=xa,
                                    rebalancing=rebalancing)
                cars.supergraph2G(tNet_cavs)
            elif alg == 'CARS3':
                if j==0:
                    cars.solve_CARS2(tNet_cavs, fcoeffs=fcoeffs, exogenous_G=False,  rebalancing=rebalancing)
                else:
                    cars.solve_CARS2(tNet_cavs, fcoeffs=fcoeffs, exogenous_G=tNet_non_cavs.G_supergraph,
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
        cavsCost.append(cavscost_/tNet_cavs.totalDemand*60)
        noncavsCost.append(noncavscost_/ tNet_non_cavs.totalDemand*60)
        totCost.append((cavscost_+noncavscost_) / (
                    tNet_cavs.totalDemand + tNet_non_cavs.totalDemand)*60)

        cavsFlow.append(cars.get_amod_flow(tNet_cavs))
        nonCavsFlow.append(tnet.get_total_G_flow(tNet_non_cavs.G))
        pedestrianFlow.append(cars.get_pedestrian_flow(tNet_cavs))
        rebalancingFlow.append(cars.get_rebalancing_flow(tNet_cavs))

        mkdir_n('results/' + dir_out)
        if rebalancing == True:
            fig0.savefig('results/' + dir_out + '/iteration-rebalance-pR-+' + str(round(penetration_rate, 1)) + '.pdf')
        else:
            fig0.savefig('results/' + dir_out + '/iteration-rebalance-pR-+' + str(round(penetration_rate, 1)) + '.pdf')

        print(penetration_rate)
        print([cavsCost, noncavsCost, totCost,])
        del tNet_cavs, tNet_non_cavs

    return cavsCost, noncavsCost, totCost, cavsFlow, nonCavsFlow, pedestrianFlow, rebalancingFlow





xa = 1.2
rebalancing = True

if rebalancing == True:
    #netFile, gFile, fcoeffs, tstamp, dir_out = tnet.get_network_parameters('EMA', experiment_name='EMA_penRate_comparison-'+'REB')
    netFile, gFile, fcoeffs, tstamp, dir_out = tnet.get_network_parameters('NYC_Uber_small_1', experiment_name='NYC_Uber_small_1_penRate_comparison-REB')

else:
    #netFile, gFile, fcoeffs, tstamp, dir_out = tnet.get_network_parameters('EMA',experiment_name='EMA_penRate_comparison-' + 'NO-REB')
    netFile, gFile, fcoeffs, tstamp, dir_out = tnet.get_network_parameters('NYC_Uber_small_1', experiment_name='NYC_Uber_small_1_penRate_comparison-NO-REB')


j = list(np.linspace(0.01, 0.99, 11))


fig1 = plt.figure(num=2, figsize=(4, 2.8))
ax1 = fig1.add_subplot()
lstyle = ['-', '--', ':']
i=0
for alg in ['disjoint']:# 'disjoint']:#, 'CARS2', 'disjoint']:

#for alg in ['disjoint']:#['CARS2','CARS']:

    cavsCost, noncavsCost, totCost, cavsFlow, nonCavsFlow, pedestrianFlow, rebalancingFlow = penRate(netFile, gFile, fcoeffs, \
                                                                                                     alg=alg, xa=xa, rebalancing=rebalancing)
    ax1.plot(j, totCost, label=alg+' (Total)', color='green', linestyle=lstyle[i], linewidth=2, marker='o')
    ax1.plot(j, cavsCost, label=alg+' (AMoDs)', color='blue', linestyle=lstyle[i], linewidth=2, marker="^")
    ax1.plot(j, noncavsCost, label=alg+' (Private)', color='red', linestyle=lstyle[i], linewidth=2, marker='x')
    ax1.legend()
    i +=1
plt.xlabel('Penetration Rate')
plt.ylabel('Avg. Travel Time (min)')
ax1.grid(True)
plt.tight_layout()
plt.xlim((0,1))
plt.ylim((18,100))
mkdir_n('results/' + dir_out)
if rebalancing ==True:
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
plt.xticks(x, x_name)
#plt.legend((p1[0], p2[0], p3[0]), ('Private Vehicle', 'AMoD', 'Rebalancing'), fontsize=12)#, 'Pedestrian'))
#plt.legend((p1[0], p2[0], p3[0], p4[0]), ('Private Vehicle', 'AMoD', 'Rebalancing', 'Micromobility'), fontsize=12)
plt.legend((p1[0], p2[0], p3[0], p4[0]), ('Private Vehicle', 'AMoD', 'Rebalancing', 'Pedestrian'), fontsize=12)
ax2.tick_params(axis='both', which='major', labelsize=11)

#plt.show()
plt.tight_layout()
mkdir_n('results/' + dir_out)

if rebalancing ==True:
    plt.savefig('results/' + dir_out +'/flows-rebalance.pdf')
else:
    plt.savefig('results/' + dir_out +'/flows-no-rebalance.pdf')

plt.show()

