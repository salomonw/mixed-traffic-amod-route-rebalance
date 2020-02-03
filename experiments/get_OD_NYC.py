

from src.utils import *
import src.tnet as tnet
odFile = "../tnet/results/joint/NYC_case_final_v2/output/NYC_OD_demand_j.csv"
g = csv2dict(odFile)
g_new = {}
for k,v in g.items():
    i,j = k.split(",")
    i = int(i.split("(")[1])
    j = int(j.split(")")[0])
    g_new[(i,j)] = float(v)
print(g_new)
tnet.writeTripsfile(g_new, 'data/trips/NYC_small_trips.txt')