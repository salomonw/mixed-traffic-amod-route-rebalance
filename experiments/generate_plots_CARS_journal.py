import src.tnet as tnet
import src.CARS as cars
from matplotlib import rc
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import pandas as pd
import experiments.build_NYC_subway_net as nyc
import random as random


rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})


g = [1,3,6]
net = 'NYC'

def genearte_plots(dir, net, g):
    fname = dir+


def plot_accuracy(fname):
    df = pd.read_csv(fname)