from src.utils import *
import src.tnet as tnet
import matplotlib.pyplot as plt
import src.CARS as cars

#netFile, gFile, fcoeffs = tnet.get_network_parameters('Braess1')
#netFile, gFile, fcoeffs = tnet.get_network_parameters('EMA')
netFile, gFile, fcoeffs = tnet.get_network_parameters('NYC_small')
xa_vec = [0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3]
#xa_vec = [0.3,0.8, 1.3]
demand_multiplier = list(np.linspace(0.6,2.5,20))


# Solve the problem for the original BPR function
real_obj = []
for g_multi in demand_multiplier:
    tNet = tnet.tNet(netFile=netFile, gFile=gFile, fcoeffs=fcoeffs)
    g_per = tnet.perturbDemandConstant(tNet.g, g_multi)
    tNet.set_g(g_per)

    tNet.solveMSAsocial()

    print(g_multi)
    socialObj = tnet.get_totalTravelTime(tNet.G, fcoeffs)
    real_obj.append(socialObj)



# Solve the problem for 3-line approx function
CARS2_obj = []
CARS2_diff = []
i=0
for g_multi in demand_multiplier:
    tNet = tnet.tNet(netFile=netFile, gFile=gFile, fcoeffs=fcoeffs)
    g_per = tnet.perturbDemandConstant(tNet.g, g_multi)
    tNet.set_g(g_per)

    tNet.build_supergraph()
    pedestrian = [(u, v) for (u, v, d) in tNet.G_supergraph.edges(data=True) if d['type'] == 'p']
    connector = [(u, v) for (u, v, d) in tNet.G_supergraph.edges(data=True) if d['type'] == 'f']

    cars.solve_CARS2(tNet, fcoeffs=fcoeffs, exogenous_G=False, rebalancing=False)
    print(g_multi)

    CARS2obj = tnet.get_totalTravelTime(tNet.G_supergraph, fcoeffs)
    CARS2_obj.append(CARS2obj)
    CARS2_diff.append((CARS2obj - real_obj[i])/real_obj[i]*100)
    i+=1

# Solve the problem for different BPR approximations and for different demand multipliers
cars_obj = {}
obj_diff = {}
flow_diff = {}
for xa in xa_vec:
    cars_obj[xa] =[]
    obj_diff[xa] = []
    flow_diff[xa] = []
    print(xa)
    i=0
    for g_multi in demand_multiplier:
        tNet = tnet.tNet(netFile=netFile, gFile=gFile, fcoeffs=fcoeffs)
        g_per = tnet.perturbDemandConstant(tNet.g, g_multi)
        tNet.set_g(g_per)

        tNet.build_supergraph()
        pedestrian = [(u, v) for (u, v, d) in tNet.G_supergraph.edges(data=True) if d['type'] == 'p']
        connector = [(u, v) for (u, v, d) in tNet.G_supergraph.edges(data=True) if d['type'] == 'f']

        cars.solve_CARS(tNet, exogenous_G=False, fcoeffs=fcoeffs, rebalancing=False, xa=xa)
        #cars.solve_CARS_noRebalancing(tNet, fcoeffs=fcoeffs, exogenous_G=False, xa=xa)
        carsObj = tnet.get_totalTravelTime(tNet.G_supergraph, fcoeffs)
        cars_obj[xa].append(carsObj)

        obj_diff[xa].append((carsObj - real_obj[i])/real_obj[i]*100)
        #flow_diff[xa].append(tnet.normFlowDiff(tNet.G, tNet.G_supergraph))
        i+=1



plt.figure()
plt.plot(demand_multiplier, real_obj, label="BPR", linestyle='--', color='k')
plt.plot(demand_multiplier, CARS2_obj, label="CARS2", linestyle='--', color='blue')
for xa in xa_vec:
    plt.plot(demand_multiplier, cars_obj[xa], label="xa="+str(xa))
plt.legend()
plt.ylabel('Obj')
plt.xlabel('Demand')


plt.figure()
plt.plot(demand_multiplier, CARS2_diff, label="CARS2", linestyle='--', color='blue')
for xa in xa_vec:
    plt.plot(demand_multiplier, obj_diff[xa], label="xa="+str(xa))
plt.legend()
plt.ylabel('Obj diff (%)')
plt.xlabel('Demand')

plt.show()
plt.figure()
for xa in xa_vec:
    plt.plot(demand_multiplier, flow_diff[xa], label="xa="+str(xa))
plt.legend()
plt.ylabel('Flow norm')
plt.xlabel('Demand')
plt.show()


plt.figure()
plt.plot(demand_multiplier, flow_diff, label='Flow norm')
plt.legend()
plt.ylabel('Obj diff')
plt.xlabel('Demand')

plt.show()
