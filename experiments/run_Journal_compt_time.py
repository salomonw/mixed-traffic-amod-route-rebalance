import src.tnet as tnet
import src.CARS as cars
import pandas as pd
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#rc('text', usetex=True)

net_names = ['NYC_Uber_small_1', 'EMA_mid', 'Anaheim']#, 'Winnipeg', 'Barcelona']#,  'Sydeny' ]
net_names = ['Anaheim']#, 'Winnipeg', 'Barcelona']#,  'Sydeny' ]

n = [2+i for i in range(6)]

#print("\ntestCars progressBar:")
#progBar = progressBar(len(n)*2*len(demand_multiplier))
#progBar.set()

v = []
for net_name in net_names:
    netFile, gFile, fcoeffs, tstamp, dir_out = tnet.get_network_parameters(net_name, experiment_name=net_name + '_comptTime_')
    tNet = tnet.tNet(netFile=netFile, gFile=gFile, fcoeffs=fcoeffs)
    tNet.solveMSA(exogenous_G=False)
    TAP_obj = tnet.get_totalTravelTime(tNet.G, fcoeffs=fcoeffs, G_exogenous=False)
    #print(TAP_obj)
    #'''
    #Table[net_name] = {}
    for i in n:
        d = {'Net':net_name, 'A':tNet.nLinks, 'V':tNet.nNodes,'W':tNet.nOD, 'n':i}
        for linear in [False, True]:
            tNet = tnet.tNet(netFile=netFile, gFile=gFile, fcoeffs=fcoeffs)
            tNet.build_supergraph(walk_multiplier=0.001)
            tNet, runtime, od_flows = cars.solve_bush_CARSn(tnet=tNet,
                                                            fcoeffs=tNet.fcoeffs,
                                                            n=i,
                                                            linear=linear,
                                                            bush=True,
                                                            rebalancing=False)
            CARSobj = tnet.get_totalTravelTime(tNet.G_supergraph, fcoeffs)
            if linear == True:
                d['t (LP)'] = runtime
                d['obj (LP)'] = (TAP_obj/CARSobj-1)*100
            else:
                d['t (QP)'] = runtime
                d['obj (QP)'] = (TAP_obj/CARSobj-1)*100
        v.append(d)
        df = pd.DataFrame(v)
        print('--------')
        print(df)
            #progBar.tic()
    del tNet
#'''

'''
fig, ax = plt.subplots(figsize=(5,2))
ax.plot(n, [v[True] for k,v in CARS.items()], label = 'LP')
ax.plot(n, [v[False] for k,v in CARS.items()], label = 'QP')
ax.set_xlabel('n')
ax.set_ylabel('% deviation from NLP')
ax.set_xlim([n[0], n[-1]])
ax.legend(framealpha=1)
ax.grid(True)
#plt.tight_layout()
plt.show()
'''



'''
print('---- solving NLP problem to set up a base ---')
real_obj = []
for g_multi in demand_multiplier:
    tNet = tnet.tNet(netFile=netFile, gFile=gFile, fcoeffs=fcoeffs)
    tNet.build_supergraph(walk_multiplier=1)
    pedestrian = [(u, v) for (u, v, d) in tNet.G_supergraph.edges(data=True) if d['type'] == 'p']
    connector = [(u, v) for (u, v, d) in tNet.G_supergraph.edges(data=True) if d['type'] == 'f']
    g_per = tnet.perturbDemandConstant(tNet.g, g_multi)
    tNet.set_g(g_per)
    cars.solve_social_Julia(tNet, exogenous_G=False)
    print('\t solve for g_multiplier = ' + str(round(g_multi,2)))
    socialObj = tnet.get_totalTravelTime(tNet.G_supergraph, fcoeffs)
    real_obj.append(socialObj)
    print(socialObj)
'''


