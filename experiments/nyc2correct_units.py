import src.tnet as tnet
import os

def change_units(G):
    for i, j in G.edges():
        G[i][j]['length'] = G[i][j]['length'] / 1600
        G[i][j]['t_0'] = G[i][j]['t_0'] / 60 / 60

def run(dir):
    gFile = "data/trips/NYC_M_trips.txt"
    fcoeffs = [1, 0, 0, 0, 0.15, 0]
    for filename in os.listdir(dir):
        netFile = os.path.join(dir, filename)
        G = tnet.readNetFile(netFile)
        #if filename.endswith('Road.txt'):
        #node_id_map = {k: G.nodes[k]['node name'] for k in G.nodes()}
        #g = tnet.readODdemandFile(gFile,node_id_map)
        #tNet = tnet.tNet(netFile=netFile, gFile=None, fcoeffs=fcoeffs)
        change_units(G)
        tnet.writeNetfile2(G, dir+'nu'+filename)

run('data/net/NYC_M/')