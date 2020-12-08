import src.tnet as tnet
import matplotlib.pyplot as plt

net = 'EMA_mid'
netFile, gFile, fcoeffs, tstamp, dir_out = tnet.get_network_parameters(net_name=net,

                                                                       experiment_name=net + 'topo_plot')
print(netFile)
tNet = tnet.tNet(netFile=netFile, gFile=gFile, fcoeffs=fcoeffs)
tNet.read_node_coordinates('data/pos/'+net+'.txt')
tnet.plot_network(tNet.G, width=0.3)
plt.show()