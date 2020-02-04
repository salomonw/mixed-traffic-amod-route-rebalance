
from src.utils import *
import src.tnet as tnet
import matplotlib.pyplot as plt
import src.CARS as cars

netFile = "data/net/EMA_net.txt"
gFile = "data/trips/EMA_trips.txt"
fcoeffs = [1,0,0,0,0.15,0]
demand_multiplier = list(np.linspace(0.5,3,25))
xa = 1


cars_obj = []
real_obj = []
diff = []
flow_diff = []
for g_multi in demand_multiplier:
    tNet = tnet.tNet(netFile=netFile, gFile=gFile, fcoeffs=fcoeffs)
    g_per = tnet.perturbDemandConstant(tNet.g, g_multi)
    tNet.set_g(g_per)

    tNet.build_supergraph()
    pedestrian = [(u, v) for (u, v, d) in tNet.G_supergraph.edges(data=True) if d['type'] == 'p']
    connector = [(u, v) for (u, v, d) in tNet.G_supergraph.edges(data=True) if d['type'] == 'f']

    cars.solve_CARS2_noRebalancing(tNet, fcoeffs=fcoeffs, exogenous_G=False, xa=xa)
    tNet.solveMSAsocial()

    print(g_multi)
    carsObj = tnet.get_totalTravelTime(tNet.G_supergraph, fcoeffs)
    socialObj = tnet.get_totalTravelTime(tNet.G, fcoeffs)
    cars_obj.append(carsObj)
    real_obj.append(socialObj)
    diff.append(abs(carsObj-socialObj))
    flow_diff.append(tnet.normFlowDiff(tNet.G, tNet.G_supergraph))

plt.figure()
plt.plot(demand_multiplier, real_obj, label="BPR")
plt.plot(demand_multiplier, cars_obj, label="Approx")
plt.plot(demand_multiplier, diff, label='Approximation error')
plt.legend()

plt.figure()
plt.plot(demand_multiplier[0:15], real_obj[0:15], label="BPR")
plt.plot(demand_multiplier[0:15], cars_obj[0:15], label="Approx")
plt.plot(demand_multiplier[0:15], diff[0:15], label='Approximation error')
plt.legend()
plt.show()


plt.figure()
plt.plot(demand_multiplier, flow_diff, label='Flow norm')
plt.legend()
plt.show()
