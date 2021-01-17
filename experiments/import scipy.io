import scipy.io
import networkx as nx

d =  scipy.io.loadmat('DataNYC.mat')
coord = d['Full'][0][0][1]
A = d['Full'][0][0][2]
