import src.tnet as tnet
import src.CARS as cars
import os
from datetime import datetime


def build_NYC_net(dir, only_road=False, symbls=False):
    gFile = "data/trips/NYC_trips.txt"
    fcoeffs = [1, 0, 0, 0, 0.15, 0]
    for filename in os.listdir(dir):
        if filename.endswith('Road_net.txt'):
            netFile = os.path.join(dir, filename)
            tNet = tnet.tNet(netFile=netFile, gFile=gFile, fcoeffs=fcoeffs)

    symbs =[]
    if only_road == False:
        tNet.G_supergraph = tNet.G.copy()
        for filename in os.listdir(dir):
            if filename.endswith(".txt") and not filename.endswith('Road_net.txt'):
                netFile = os.path.join(dir, filename)
                layer = tnet.readNetFile(netFile=netFile)
                if 'Walk' in filename:
                    layer_symb="'"
                else:
                    layer_symb=filename.split('_')[3]
                    symbs.append(layer_symb)
                tNet.add_layer(layer=layer, layer_symb=layer_symb)
            else:
                continue



    tstamp = datetime.now().strftime("%Y-%m-%d_%H:%M:%S")
    fcoeffs = [1, 0, 0, 0, 0.15, 0]
    if symbls == True:
        return tNet, tstamp, fcoeffs, symbs
    else:
        return tNet, tstamp, fcoeffs




