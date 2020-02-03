import networkx as nx
from gurobipy import *
from src.utils import *
import numpy as np
import src.msa as msa
import matplotlib.pyplot as plt
import src.trafficAssignment.assign as ta
from scipy.special import comb
from scipy.linalg import block_diag
import copy


class tNet():
    def __init__(self, netFile=None, gFile=None, G=None, g=None, fcoeffs=[1, 0, 0, 0, 0.15, 0]):
        """
        Initialize a Traffic Network object. It requieres networkx and Gurobi libraries

        Parameters
        ----------
        self : a tNet object

        netFile: a network file in the format proposed by Bar Gera

        gFile: an OD demand file in the format proposed by Bar Gera

        Returns
        -------
        A tNet Object.

        """

        if G == None and netFile != None:
            G = readNetFile(netFile)
        elif G != None and netFile == None:
            G = G

        node_id_map = {k: G.nodes[k]['node name'] for k in G.nodes()}
        if g == None and gFile != None:
            g = readODdemandFile(gFile, node_id_map)
        elif g != None and gFile == None:
            g = g

        self.netFileName = netFile
        self.gFileName = gFile
        self.nLinks = len(G.edges())
        self.nNodes = len(G.nodes())
        self.nOD = len(g)
        self.Zones = getZones(g)
        self.nZones = getNumZones(g)
        self.totalDemand = sum(g.values())
        self.g = g
        self.gGraph = buildDemandGraph(g)
        self.G = G
        self.node_id_map = node_id_map
        self.fcoeffs = fcoeffs
        self.nPoly = len(fcoeffs)
        self.TAP = self.build_TAP()
        self.incidence_matrix, self.link_id_dict = incidence_matrix(self.G)

    def build_TAP(self):
        """
	    Build a traffic assignment object based on the traffic network 
	    Jurgen Hackl <hackl@ibi.baug.ethz.ch>

	    Parameters
	    ----------

		gdict: OD demands dict

	    Returns
	    -------
	    An nx object.

	    """
        assert (self.nZones > 0), "Number of Zones is zero!!"
        assert (self.totalDemand > 0), "Total demand is zero!!"
        assert (self.nNodes > 0), "No nodes in graph!!"
        assert (self.nLinks > 0), "No links in graph!!"
        TAP = msa.TrafficAssignment(self.G, self.gGraph, fcoeffs=self.fcoeffs, iterations=350)
        return TAP

    def build_pedestrian_net(self):
        """
        build a pedestrian net by considering both directions links

        Parameters
        ----------

        self: a tnet object

        Returns
        -------
        An attribute on the object containing the pedestrian network

        """
        G = nx.DiGraph()
        for i, j in self.G.edges():
            G.add_edge(i, j, length=self.G[i][j]['length'])
            G.add_edge(j, i, length=self.G[i][j]['length'])
        self.G_pedestrian = G

    def build_supergraph(self, walk_multiplier=7):
        """
        build a supergraph mixing pedestrian and vehicle networks

        Parameters
        ----------

        self: a tnet object

        Returns
        -------
        An attribute on the object containing the pedestrian network

        """
        G = self.G.copy()
        for i, j in self.G.edges():
            G.add_edge(str(i) + "'", str(j) + "'", length=self.G[i][j]['length'],
                       t_0=self.G[i][j]['length'] * walk_multiplier, capacity=10000000000, type='p')
            G.add_edge(str(j) + "'", str(i) + "'", length=self.G[i][j]['length'],
                       t_0=self.G[i][j]['length'] * walk_multiplier, capacity=10000000000, type='p')
            G.add_edge(i, str(i) + "'", t_0=0, capacity=99999999, type='f')
            G.add_edge(str(i) + "'", i, t_0=0, capacity=99999999, type='f')
            G.add_edge(j, str(j) + "'", t_0=0, capacity=99999999, type='f')
            G.add_edge(str(j) + "'", j, t_0=0, capacity=99999999, type='f')
        self.G_supergraph = G

    def solveMSA(self, exogenous_G=False):
        """
	    Solve the MSA flows for a traffic network using the MSA module by 
	    Jurgen Hackl <hackl@ibi.baug.ethz.ch>

	    Parameters
	    ----------

		gdict: OD demands dict

	    Returns
	    -------
	    An nx object.

	    """
        self.TAP.run(fcoeffs=self.fcoeffs, build_t0=False, exogenous_G=exogenous_G)
        self.G = self.TAP.graph

    def solveMSAsocial(self, build_t0=False, exogenous_G=False):
        """
	    Solve the MSA social flows for a traffic network using the MSA module by
	    Jurgen Hackl <hackl@ibi.baug.ethz.ch>

	    Parameters
	    ----------

		gdict: OD demands dict

	    Returns
	    -------
	    An nx object.

	    """
        self.TAP.run_social(fcoeffs=self.fcoeffs, build_t0=build_t0, exogenous_G=exogenous_G)
        self.G = self.TAP.graph

    def set_fcoeffs(self, fcoeff):
        """
	    set function coefficients of the transportation netwrok object

	    Parameters
	    ----------

		fcoeffs: function coefficients of a polynomial function

	    Returns
	    -------
	    updates the fcoeffs attribute of the tnet object

	    """
        self.fcoeffs = fcoeff

    # TAP = msa.TrafficAssignment(self.G, self.gGraph, fcoeffs=fcoeff)

    def set_g(self, g):
        """
	    sets a demand dictionary and updates the the tNet and TAP objects

	    Parameters
	    ----------

		g: an OD demand dictionary

	    Returns
	    -------
	    updates the objects in the tnet where g takes place

	    """
        self.nOD = len(g)
        self.Zones = getZones(g)
        self.nZones = getNumZones(g)
        self.totalDemand = sum(g.values())
        self.g = g
        self.gGraph = buildDemandGraph(g)
        self.TAP = self.build_TAP()

    def read_flow_file(self, fname):
        """
	    Read flow file in tntp format. If G is provided, 
	    we add the info to G, otherwise, generate a nx object.

	    Parameters
	    ----------
		
		fname: input file name
	    G: networkx obj

	    Returns
	    -------
	    a networkx object 

	    """
        read_flowFile(fname, self.G)

    def write_flow_file(self, fname):
        """
	    Write flow file in tntp format. If G is provided, 
	    we add the info to G, otherwise, generate a nx object.

	    Parameters
	    ----------
		
		fname: input file name
	    G: networkx obj

	    Returns
	    -------
	    a networkx object 

	    """
        write_flowFile(self.G, fname)

    def travel_time(self, i, j):
        """
        evalute the travel time function for edge i->j

        Parameters
        ----------

        i: starting node of edge
        j: ending node of edge

        Returns
        -------
        float

        """
        return sum([self.fcoeffs[n] * (self.G[i][j]['flow'] / self.G[i][j]['capacity']) ** n for n in
                    range(len(self.fcoeffs))])

    def eval_obj(self):
        """
        evalute total travel time in network

        Parameters
        ----------
        self

        Returns
        -------
        float

        """
        return sum([self.G[i][j]['flow'] * self.G[i][j]['t_0'] * self.travel_time(i, j) for i, j in self.G.edges()])

    def read_node_coordinates(self, fname):
        """
        reads node coordinates to the graph network

        Parameters
        ----------
        self
        fname: file containing the node coordinates

        Returns
        -------
        modified object

        """
        with open(fname) as file:
            file_lines = file.readlines()

        od_d = {}
        od_d_t = {}
        for line in file_lines:
            if "x" in line:
                continue
            if ";" in line:
                line_ = line.split()
                n = int(line_[0])
                x = float(line_[1])
                y = float(line_[2])
                self.G.nodes[n]['pos'] = (x, y)


def plot_network_flows(G, weight='flow', width=3, cmap=plt.cm.Blues):
    # TODO: add explaination
    fig, ax = plt.subplots()
    pos = nx.get_node_attributes(G, 'pos')
    edges, weights = zip(*nx.get_edge_attributes(G, weight).items())
    nx.draw(G, pos, node_color='b', edgelist=edges, edge_color=weights, width=width, edge_cmap=cmap)
    labels = {(i, j): int(G[i][j][weight]) for i, j in G.edges()}
    nx.draw_networkx_edge_labels(G, pos=pos, edge_labels=labels)
    return fig, ax


def perturbDemandConstant(g, constant):
    """
    Perturb demand by a random constant value

    Parameters
    ----------

	g: demand dict
	max_var: maximum percentage of the demands which is goinf to variate either
			increase or decrease

    Returns
    -------
    a perturbed demand dict

    """
    a = {}
    for od, demand in g.items():
        a[od] = demand * constant
    return a


def incidence_matrix(G):
    """
    build incidence matrix and column index dictionary. Note that since 
    node 0 don't exist, the node is the column number plus one

    Parameters
    ----------

	a nx element

    Returns
    -------
    A sparse matrix and a dictionary to match line number and link

    """
    nNodes = len(G.nodes())
    nLinks = len(G.edges())
    N = np.zeros((nNodes, nLinks))
    link_dict = {}
    idx = 0
    for s, t in G.edges():
        link_dict[idx] = (s, t)
        N[s - 1, idx] = -1
        N[t - 1, idx] = 1
        idx += 1
    return N, link_dict


def buildDemandGraph(g):
    """
    a nx graph defininf OD demands

    Parameters
    ----------

	gdict: OD demands dict

    Returns
    -------
    An nx object.

    """
    od_graph = nx.DiGraph()
    for (s, t), d in g.items():
        od_graph.add_edge(s, t, demand=d)
    return od_graph


def readNetFile(netFile, sep="\t"):
    """
    Read the netFile and convert it to a nx object

    Parameters
    ----------

    netFile: a network file in the format proposed by Bar Gera

    Returns
    -------
    An nx object.

    """

    # Create a networkx obj
    G = nx.DiGraph()
    # Read the network file
    with open(netFile) as file_flow:
        file_flow_lines = file_flow.readlines()

    for line in file_flow_lines:
        if ";" in line and "~" not in line:
            links = line.split(sep)
            G.add_edge(int(links[1]), int(links[2]), capacity=float(links[3]), \
                       length=float(links[4]), t_0=float(links[5]), \
                       B=float(links[6]), power=float(links[7]), speedLimit=float(links[8]), \
                       toll=float(links[9]), type=float(links[10]))
    G = nx.convert_node_labels_to_integers(G, first_label=1, ordering='sorted', label_attribute='node name')

    return G


def readODdemandFile(gfile, node_id_map):
    """
	Read the gfile and convert it to a dict

	Parameters
	----------

	gFile: a demand file in the format proposed by Bar Gera

	Returns
	-------
	A dict with OD as key and demand as value

	"""

    with open(gfile) as trips:
        trip_lines = trips.readlines()

    od_d = {}
    od_d_t = {}
    for line in trip_lines:
        if "Origin" in line:
            origin = node_id_map[int(line.split("gin")[1])]
        if ";" in line:
            line_ = line.split(";")
            for j in line_:
                if ":" in j:
                    dest = node_id_map[int(j.split(":")[0])]
                    d = float(j.split(":")[1])
                    if origin != dest:
                        od_d[(origin, dest)] = d
    return od_d


@timeit
def writeNetfile(G, g, fname):
    """
	write net file from G (networkx)

	Parameters
	----------

	G: a neteworkx object

	Returns
	-------
	file
	"""
    nZones = str(len(set([i for i, j in g.keys()])))
    nNodes = str(len(G.nodes()))
    nLinks = str(len(G.edges()))
    header = "<NUMBER OF ZONES> " + nZones + "\n<NUMBER OF NODES> " + nNodes + "\n<FIRST THRU NODE> 1\n<NUMBER OF LINKS> " + nLinks + "\n<END OF METADATA>\n~  Init node  Term node  Capacity  Length  Free Flow Time  B  Power  Speed limit  Toll  Type  ;\n"
    text = ""
    idx = 0
    link_id = {}
    for (s, t) in G.edges():
        idx += 1
        link_id[idx] = (s, t)
        text += "\t" + str(s) + "\t" + str(t) + "\t" + str(G[s][t]["capacity"]) + "\t" + str(
            G[s][t]["length"]) + "\t" + str(G[s][t]["t_0"]) \
                + "\t" + str(G[s][t]["B"]) + "\t" + str(G[s][t]["power"]) + "\t" + str(
            G[s][t]["speedLimit"]) + "\t" + str(G[s][t]["toll"]) \
                + "\t" + str(G[s][t]["type"]) + "\t;\n"
    write_file(header + text, fname)
    return link_id, header + text


@timeit
def writeTripsfile(g, fname):
    """
	write trips file from dict 

	Parameters
	----------
]
	g dict

	Returns
	-------
	file
	"""
    nZones = str(len(set([i for i, j in g.keys()])))
    totalFlow = sum([d for d in g.values()])
    header = "<NUMBER OF ZONES> " + nZones + "\n<TOTAL OD FLOW> " + str(totalFlow) + "\n<END OF METADATA>\n\n"

    text = ""
    nodes = list(set([s for s, t in g.keys()]))
    nodes.extend(list(set([t for s, t in g.keys()])))
    nodes = list(set(nodes))
    for o in nodes:
        txt = ""
        txt = "Origin " + str(o) + "\n"
        demandTxt = ""
        for d in nodes:
            try:
                gk = str(g[(o, d)])
            except:
                gk = str(0)
            demandTxt += str(d) + "\t : \t" + str(gk) + ";\n"
        text += txt + demandTxt + "\n\n"
    write_file(header + text, fname)
    return header + text


def getZones(gdict):
    """
    Returns Zones in a OD file

    Parameters
    ----------

    gdict: a demand dictiornary

    Returns
    -------
    a list with the Zones of the network

    """
    sources = [s for (s, t) in gdict.keys()]
    targets = [t for (s, t) in gdict.keys()]
    sources.extend(targets)
    return set(sources)


def getNumZones(gdict):
    """
    Finds the number of zones in a network

    Parameters
    ----------

    gdict: a demand dictiornary

    Returns
    -------
    number of zones

    """
    return len(getZones(gdict))


def greenshieldFlow(speed, capacity, free_flow_speed):
    """
    Returns the flow of a link following the Fundamental Diagram (Greenshield's Model).

    Parameters
    ----------

    speed: miles/hr or km/hr
    capacity: link's capacity in veh/hr
    free_flow_speed: link's free flow speed in miles/hr or km/hr

    Returns
    -------
    flow: resulting flow 

    """
    if speed > free_flow_speed or capacity < 0:
        return 0
    x = 4 * capacity * speed / free_flow_speed - 4 * capacity * (speed ** 2) / (free_flow_speed ** 2)
    return x


def travel_time(G, fcoeffs, i, j, G_exogenous=False):
    """
	evalute the travel time function for edge i->j

	Parameters
	----------
	tnet: transportation network object
	i: starting node of edge
	j: ending node of edge

	Returns
	-------
	float

	"""
    if G_exogenous:
        return sum([fcoeffs[n] * ((G[i][j]['flow'] + G_exogenous[i][j]['flow']) / G[i][j]['capacity']) ** n for n in
                    range(len(fcoeffs))])
    else:
        return sum([fcoeffs[n] * (G[i][j]['flow'] / G[i][j]['capacity']) ** n for n in range(len(fcoeffs))])


def get_totalTravelTime(G, fcoeffs, G_exogenous=False):
    """
    Return the total travel time of a network

    Parameters
    ----------

    G: networkx obj

    Returns
    -------
    a float 

    """
    return sum([G[i][j]['flow'] * G[i][j]['t_0'] * travel_time(G, fcoeffs, i, j, G_exogenous=G_exogenous) for i, j in
                G.edges()])


def read_flowFile(fname, G=False):
    """
	Read flow file in tntp format. If G is provided, 
	we add the info to G, otherwise, generate a nx object.

	Parameters
	----------

	fname: input file name
	G: networkx obj

	Returns
	-------
	a networkx object 

	"""
    with open(fname) as flows:
        flow_lines = flows.readlines()
    if G == False:
        G = nx.DiGraph()
    for line in flow_lines:
        if "Tail" not in line and ";" in line:
            line_ = line.split("\t")
            i = 0
            for j in line_:
                if i == 1:
                    s = int(j)
                elif i == 2:
                    t = int(j)
                elif i == 3:
                    flow = float(j)
                i += 1
            if G == False:
                G.add_edge(s, t)
                G[s][t]['flow'] = flow
            else:
                G[s][t]['flow'] = flow
    return G


def write_flowFile(G, fname):
    """
    Write flow file in tntp format. If G is provided, 
    we add the info to G, otherwise, generate a nx object.

    Parameters
    ----------
	
	fname: input file name
    G: networkx obj

    Returns
    -------
    a networkx object 

    """
    nNodes = str(len(G.nodes()))
    nLinks = str(len(G.edges()))
    header = "~\tTail\tHead\tVolume\t;\n"
    text = ""
    idx = 0
    link_id = {}
    for (s, t) in G.edges():
        idx += 1
        link_id[idx] = (s, t)
        text += "\t" + str(s) + "\t" + str(t) + "\t" + str(G[s][t]["flow"]) + "\t;\n"
    write_file(header + text, fname)


def solveMSA_julia(tnet, fcoeffs=False, G_exo=False, net_fname="tmp_jl/net.txt", trips_fname="tmp_jl/trips.txt"):
    """
    Solve the MSA flows for a traffic network using the MSA module by
    __ .. __ in Julia

    Parameters
    ----------

	tnet object

    Returns
    -------
    An updated nx object.
    """
    if fcoeffs == False:
        fcoeffs = tnet.fcoeffs
    # pwd = os.getcwd()
    # link_id, text_net = writeNetfile(tnet.G, tnet.g, net_fname)
    # text_trips = writeTripsfile(tnet.g, trips_fname)
    new_G = ta.assignment(tnet.G, tnet.g, tnet.fcoeffs, G_exo=G_exo, flow=False, method='MSA', accuracy=0.0001,
                          max_iter=2000)
    tnet.G = new_G

    return tnet


def solveMSA_julia_social(tnet, fcoeffs=False, G_exo=False, net_fname="tmp_jl/net.txt", trips_fname="tmp_jl/trips.txt"):
    """
    Solve the MSA flows for a traffic network using the MSA module by
    __ .. __ in Julia

    Parameters
    ----------

	tnet object

    Returns
    -------
    An updated nx object.
    """
    if fcoeffs == False:
        fcoeffs = tnet.fcoeffs
    # pwd = os.getcwd()
    # link_id, text_net = writeNetfile(tnet.G, tnet.g, net_fname)
    # text_trips = writeTripsfile(tnet.g, trips_fname)
    new_G = ta.assignment_social(tnet.G, tnet.g, tnet.fcoeffs, G_exo=G_exo, flow=False, method='MSA', accuracy=0.0001,
                                 max_iter=000)
    tnet.G = new_G

    return tnet


def add_net_flows(array):
    # TODO: add description

    G = array[0].G.copy()
    for tn in array[1:]:
        for i, j in G.edges():
            G[i][j]['flow'] += tn.G[i][j]['flow']
    return G


def add_G_flows(array):
    # TODO: add description

    G = array[0].copy()
    for tn in array[1:]:
        for i, j in G.edges():
            G[i][j]['flow'] += tn[i][j]['flow']
    return G


def flowDiff(graph1, graph2):
    """
    Return the difference between the flows of two graphs.

    Parameters
    ----------

    graph1: netwrorkx obj
    graph2: networkx obj

    Returns
    -------
    a dictionary with name od edge and diff in flows. graph1-graph2

    """
    flowDiff_ = {}
    for edge in graph1.edges():
        flow1 = graph1[edge[0]][edge[1]]['flow']
        flow2 = graph2[edge[0]][edge[1]]['flow']
        flowDiff_[edge] = flow1-flow2
    return flowDiff_