from src.utils import *
import src.tnet as tnet
import src.CARS as cars
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)


#netFile, gFile, fcoeffs, tstamp, dir_out = tnet.get_network_parameters('Braess1', experiment_name='Braess1_penRate_disjoint')
#netFile, gFile, fcoeffs, tstamp, dir_out = tnet.get_network_parameters('NYC', experiment_name='NYC_CARS_precision')
#netFile, gFile, fcoeffs, tstamp, dir_out = tnet.get_network_parameters('EMA', experiment_name='EMA_CARS_precision')
netFile, gFile, fcoeffs, tstamp, dir_out = tnet.get_network_parameters('NYC_Uber_small_1', experiment_name='NYC_Uber_small_1_precision')

xa_vec = [0.6, 0.8, 1.0, 1.2]
#xa_vec = [0.6, 0.8, 1.3]
demand_multiplier = list(np.linspace(0.5,1.5,11))


# Solve the problem for the original BPR function
#'''
print('---- solving NLP problem to set up a base ---')
real_obj = []
for g_multi in demand_multiplier:
    tNet = tnet.tNet(netFile=netFile, gFile=gFile, fcoeffs=fcoeffs)
    tNet.build_supergraph()
    pedestrian = [(u, v) for (u, v, d) in tNet.G_supergraph.edges(data=True) if d['type'] == 'p']
    connector = [(u, v) for (u, v, d) in tNet.G_supergraph.edges(data=True) if d['type'] == 'f']
    g_per = tnet.perturbDemandConstant(tNet.g, g_multi)
    tNet.set_g(g_per)

    cars.solve_social_Julia(tNet, exogenous_G=False)

    print('\t solve for g_multiplier = ' + str(round(g_multi,2)))
    socialObj = tnet.get_totalTravelTime(tNet.G_supergraph, fcoeffs)
    real_obj.append(socialObj)
#'''
#real_obj = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]


# Solve the problem for 3-line approx function
print('---- solving problem with CARS3 ---')
CARS2_obj = []
CARS2_diff = []
i=0
for g_multi in demand_multiplier:
    tNet = tnet.tNet(netFile=netFile, gFile=gFile, fcoeffs=fcoeffs)
    tNet.build_supergraph()
    pedestrian = [(u, v) for (u, v, d) in tNet.G_supergraph.edges(data=True) if d['type'] == 'p']
    connector = [(u, v) for (u, v, d) in tNet.G_supergraph.edges(data=True) if d['type'] == 'f']
    g_per = tnet.perturbDemandConstant(tNet.g, g_multi)
    tNet.set_g(g_per)

    #cars.solve_CARS2(tNet, fcoeffs=fcoeffs, exogenous_G=False, rebalancing=False)
    cars.solve_bush_CARSn(tNet, fcoeffs, n=3, rebalancing=False, bush=True)
    print('\t solve for g_multiplier = ' + str(round(g_multi,2)))

    CARS2obj = tnet.get_totalTravelTime(tNet.G_supergraph, fcoeffs)
    CARS2_obj.append(CARS2obj)
    CARS2_diff.append((CARS2obj - real_obj[i])/real_obj[i]*100)
    i+=1

# Solve the problem for different BPR approximations and for different demand multipliers
print('---- solving problem with CARS2 ---')
cars_obj = {}
obj_diff = {}
flow_diff = {}
for xa in xa_vec:
    print('\t -- For theta = '+ str(xa))
    cars_obj[xa] =[]
    obj_diff[xa] = []
    flow_diff[xa] = []
    i=0
    for g_multi in demand_multiplier:
        tNet = tnet.tNet(netFile=netFile, gFile=gFile, fcoeffs=fcoeffs)
        tNet.build_supergraph()
        pedestrian = [(u, v) for (u, v, d) in tNet.G_supergraph.edges(data=True) if d['type'] == 'p']
        connector = [(u, v) for (u, v, d) in tNet.G_supergraph.edges(data=True) if d['type'] == 'f']
        g_per = tnet.perturbDemandConstant(tNet.g, g_multi)
        tNet.set_g(g_per)

        cars.solve_CARS(tNet, fcoeffs, rebalancing=False, xa=xa)

        carsObj = tnet.get_totalTravelTime(tNet.G_supergraph, fcoeffs)
        cars_obj[xa].append(carsObj)

        obj_diff[xa].append((carsObj - real_obj[i])/real_obj[i]*100)
        #flow_diff[xa].append(tnet.normFlowDiff(tNet.G, tNet.G_supergraph))
        i+=1
        print('\t \t solve for g_multiplier = ' + str(round(g_multi, 2)))


markers = ['>', '^', 'v', '<', 'o', 'x', ]
plt.figure(num=None, figsize=(6, 2.25))
#plt.plot(demand_multiplier, real_obj, label="BPR", linestyle='--', color='k')
plt.plot(demand_multiplier, CARS2_obj, label="CARS3", linestyle=':', color='black', linewidth=2)
i=0
for xa in xa_vec:
    plt.plot(demand_multiplier, cars_obj[xa], label="CARS, $\\theta=$"+str(xa), linewidth=2, marker=markers[i])
    i +=1
plt.legend()
plt.ylabel('Obj')
plt.xlabel('Demand mulitplier')
plt.tight_layout()
plt.grid(True)
plt.ylim((0,7))
plt.xlim((0.5,1.5))



plt.figure(num=None, figsize=(6, 2.25))
plt.plot(demand_multiplier, CARS2_diff, label="CARS3", linestyle=':', color='black', linewidth=2)
i=0
for xa in xa_vec:
    plt.plot(demand_multiplier, obj_diff[xa], label="CARS, $\\theta=$"+str(xa), linewidth=2, marker=markers[i])
    i+=1
plt.legend()
plt.ylabel('Obj diff (\\%)')
plt.xlabel('Demand mulitplier')
plt.grid(True)
plt.ylim((0,7))
plt.xlim((0.5,1.5))
plt.tight_layout()



mkdir_n('results/' + dir_out)
plt.savefig('results/' + dir_out +'/CARS3_precision.pdf')
plt.show()

'''
plt.show()
plt.figure()
for xa in xa_vec:
    plt.plot(demand_multiplier, flow_diff[xa], label="$\\theta=$"+str(xa))
plt.legend()
plt.ylabel('Flow norm')
plt.xlabel('Demand')
plt.show()


plt.figure()
plt.plot(demand_multiplier, flow_diff, label='Flow norm')
plt.legend()
plt.ylabel('Obj diff')
plt.xlabel('Demand')

mkdir_n('results/' + dir_out)
plt.savefig('results/' + dir_out +'/CARS3_precision.png', dpi=300)
plt.show()
'''

