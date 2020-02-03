
from src.utils import *
import src.tnet as tnet
import matplotlib.pyplot as plt
import src.CARS as cars

netFile = "data/net/EMA_net.txt"
gFile = "data/trips/EMA_trips.txt"
fcoeffs = [1,0,0,0,0.15,0]
demand_multiplier = list(np.linspace(0.5,3,25))


cars_obj = []
real_obj = []
diff = []
for g_multi in demand_multiplier:
    tNet = tnet.tNet(netFile=netFile, gFile=gFile, fcoeffs=fcoeffs)
    tNet.build_supergraph()
    pedestrian = [(u, v) for (u, v, d) in tNet.G_supergraph.edges(data=True) if d['type'] == 'p']
    connector = [(u, v) for (u, v, d) in tNet.G_supergraph.edges(data=True) if d['type'] == 'f']

    g_per = tnet.perturbDemandConstant(tNet.g, g_multi)
    tNet.set_g(g_per)

    cars.solve_CARS2(tNet, fcoeffs=fcoeffs, exogenous_G=False, xa=1.2)

    print(g_multi)
    cars_obj.append(cars.get_totalTravelTime_approx(tNet, fcoeffs=fcoeffs, xa=1.2))
    real_obj.append(tnet.get_totalTravelTime(tNet.G_supergraph, fcoeffs))
    diff.append(abs(tnet.get_totalTravelTime(tNet.G_supergraph, fcoeffs)-cars.get_totalTravelTime_approx(tNet, fcoeffs=fcoeffs, xa=1.2)))

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