import src.tnet as tnet
import matplotlib.pyplot as plt

def plot_topology(netname):
    netFile, gFile, fcoeffs, tstamp, dir_out = tnet.get_network_parameters(net_name=netname,

                                                                           experiment_name=netname + 'topo_plot')
    tNet = tnet.tNet(netFile=netFile, gFile=gFile, fcoeffs=fcoeffs)
    tNet.read_node_coordinates('data/pos/'+netname+'.txt')
    fig, ax = tnet.plot_network(tNet.G, width=0.3)
    return fig, ax


def plot_convergance(fname_sys, fname_usr):




def plot_