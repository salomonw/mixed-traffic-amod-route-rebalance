import src.tnet as tnet
import src.CARS as cars
import numpy as np
import matplotlib.pyplot as plt
import copy


def penRate(netFile, gFile, fcoeffs, alg='CARS2', rebalancing=True, xa=1):
    cavsCost = []
    noncavsCost = []
    totCost = []
    cavsFlow = []
    nonCavsFlow = []
    pedestrianFlow = []
    rebalancingFlow = []
    for penetration_rate in np.linspace(0.01, 0.99, 2):
        tNet_cavs = tnet.tNet(netFile=netFile, gFile=gFile, fcoeffs=fcoeffs)
        tNet_non_cavs = tnet.tNet(netFile=netFile, gFile=gFile, fcoeffs=fcoeffs)
        g_cavs = tnet.perturbDemandConstant(tNet_cavs.g, constant=penetration_rate)
        g_non_cavs = tnet.perturbDemandConstant(tNet_cavs.g, constant=(1 - penetration_rate))
        tNet_cavs.set_g(g_cavs)
        tNet_non_cavs.set_g(g_non_cavs)
        tNet_cavs.build_supergraph(walk_multiplier=999999)

        it = []
        for i in range(7):

            if i == 0:
                if alg == 'CARS':
                    cars.solve_CARS(tNet_cavs, exogenous_G=False, fcoeffs=fcoeffs, xa=xa, rebalancing=rebalancing)
                    cars.supergraph2G(tNet_cavs)
                elif alg == 'CARS2':
                    cars.solve_CARS2(tNet_cavs, exogenous_G=False, fcoeffs=fcoeffs, xa=xa, rebalancing=rebalancing)
                    cars.supergraph2G(tNet_cavs)
                elif alg == 'disjoint':
                    cars.solveMSAsocialCARS(tNet_cavs, exogenous_G=False)
                    cars.solve_rebalancing(tNet_cavs, exogenous_G=0)
                    #cars.G2supergraph(tNet_cavs)
                else:
                    print('please select an algorithm (CARS, CARS2 or disjoint) !!!!!!!!!')
            else:
                if alg == 'CARS':
                    cars.solve_CARS(tNet_cavs, exogenous_G=tNet_non_cavs.G, fcoeffs=fcoeffs, xa=xa, rebalancing=rebalancing)
                    cars.supergraph2G(tNet_cavs)
                elif alg == 'CARS2':
                    cars.solve_CARS2(tNet_cavs, exogenous_G=tNet_non_cavs.G, fcoeffs=fcoeffs, xa=xa, rebalancing=rebalancing)
                    cars.supergraph2G(tNet_cavs)
                elif alg == 'disjoint':
                    cars.solveMSAsocialCARS(tNet_cavs, exogenous_G=tNet_non_cavs.G)
                    cars.solve_rebalancing(tNet_cavs, exogenous_G=tNet_non_cavs.G)
                else:
                    print('please select an algorithm (CARS, CARS2 or disjoint) !!!!!!!!!')
                    break

            tNet_non_cavs.solveMSA(exogenous_G=tNet_cavs.G)

            totalCost = tnet.get_totalTravelTime(tNet_non_cavs.G, fcoeffs,
                                                 G_exogenous=tNet_cavs.G) + cars.get_totalTravelTime_without_Rebalancing(
                                                tNet_cavs, G_exogenous=tNet_non_cavs.G)
            it.append(totalCost)

        if penetration_rate > .4 and penetration_rate < 0.58:
            plt.figure()
            plt.plot(it)
            plt.show()

        cavsCost.append(cars.get_totalTravelTime_without_Rebalancing(tNet_cavs,
                                                                     G_exogenous=tNet_non_cavs.G) / tNet_cavs.totalDemand)
        noncavsCost.append(
            tnet.get_totalTravelTime(tNet_non_cavs.G, fcoeffs, G_exogenous=tNet_cavs.G) / tNet_non_cavs.totalDemand)
        totCost.append(totalCost / (
                    tNet_cavs.totalDemand + tNet_non_cavs.totalDemand))  # TODO: calculate the cost with rebalancing

        cavsFlow.append(cars.get_amod_flow(tNet_cavs))
        nonCavsFlow.append(tnet.get_total_G_flow(tNet_non_cavs.G))
        pedestrianFlow.append(cars.get_pedestrian_flow(tNet_cavs))
        rebalancingFlow.append(cars.get_rebalancing_flow(tNet_cavs))
        print(penetration_rate)
        del tNet_cavs, tNet_non_cavs

    return cavsCost, noncavsCost, totCost, cavsFlow, nonCavsFlow, pedestrianFlow, rebalancingFlow





#netFile, gFile, fcoeffs = tnet.get_network_parameters('Braess1')
#netFile, gFile, fcoeffs = tnet.get_network_parameters('EMA')
netFile, gFile, fcoeffs = tnet.get_network_parameters('NYC_small')
xa=1.2


j = list(np.linspace(0.01, 0.99, 2))

lstyle = ['-', '--', ':']
i=0
for alg in ['CARS2', 'CARS']:#, 'CARS2', 'disjoint']:

#for alg in ['disjoint']:#['CARS2','CARS']:

    cavsCost, noncavsCost, totCost, cavsFlow, nonCavsFlow, pedestrianFlow, rebalancingFlow = penRate(netFile, gFile, fcoeffs, alg=alg, xa=xa, rebalancing=False)

    plt.plot(j, totCost, label=alg+' (Total)', color='green', linestyle=lstyle[i])
    plt.plot(j, cavsCost, label=alg+' (AMoDs)', color='blue', linestyle=lstyle[i])
    plt.plot(j, noncavsCost, label=alg+' (Private)', color='red', linestyle=lstyle[i])
    plt.legend()
    i +=1

plt.show()


plt.figure()
width = 0.3
x = list(range(2))
p1 = plt.bar(x, nonCavsFlow, width)
p2 = plt.bar(x, cavsFlow, width,
             bottom=nonCavsFlow)
p3 = plt.bar(x, rebalancingFlow, width,
             bottom=[cavsFlow[i] + nonCavsFlow[i] for i in range(len(cavsFlow))])
p4 = plt.bar(x, pedestrianFlow, width,
             bottom=[cavsFlow[i] + nonCavsFlow[i] + rebalancingFlow[i] for i in range(len(cavsFlow))])
plt.ylabel('Flow')
plt.xticks(j)
plt.legend((p1[0], p2[0], p3[0], p4[0]), ('NonCAVs', 'CAVs', 'Rebalancing', 'Pedestrian'))

plt.show()


# TODO: check social solution and if they are solving for their own fleet