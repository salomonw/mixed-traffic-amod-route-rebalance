#!/usr/bin/python -tt
# -*- coding: utf-8 -*-
# =============================================================================
# Time-stamp: <Don 2018-10-11 11:16 juergen>
# File      : msa.py
# Creation  : 08 Oct 2015
#
# Copyright (c) 2015 Jürgen Hackl <hackl@ibi.baug.ethz.ch>
#               http://www.ibi.ethz.ch
# $Id$
#
# Description : A simple traffic model based on MSA and BPR
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# =============================================================================
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import itertools
import collections
import ast
from scipy.sparse import coo_matrix
#from joblib import Parallel, delayed
from multiprocessing import Pool
from src.utils import *


import multiprocessing as mp



class TrafficAssignment(object):
    """Traffic model to estimate the flow, velocity and travel time on a road
    network.
    Calculating the flow and other parameters for a given network using the
    'Method of Successive Averages' [1]_.
    Iterative algorithms were developed, at least partially, to overcome the
    problem of allocating too much traffic to low-capacity links. In an
    iterative assignment algorithm the ‘current’ flow on a link is calculated
    as a linear combination of the current flow on the previous iteration and
    an auxiliary flow resulting from an all-or-nothing assignment in the
    present iteration. The algorithm can be described by the following steps:
    1. Select a suitable initial set of current link costs, usually free-flow
       travel times. Initialise all flows :math:`V_a = 0`; make :math:`n = 0`.
    2. Build the set of minimum cost trees with the current costs; make
       :math:`n = n + 1`.
    3. Load the whole of the matrix :math:`T` all-or-nothing to these trees
       obtaining a set of auxiliary flows :math:`F_a`.
    4. Calculate the current flows as:
           .. math::
              V_a^n = (1 − \\phi)V_a^{n-1} + \\phi F_a
    5. Calculate a new set of current link costs based on the flows
       :math:`V_a^n` . If the flows (or current link costs) have not changed
       significantly in two consecutive iterations, stop; otherwise proceed to
       step 2. Another, less good but quite common, criterion for stopping is
       simply to fix the maximum number of iterations; should be calculated in
       this case as well to know how close the solution is to Wardrop’s
       equilibrium.
    Iterative assignment algorithms differ in the method used to give a value
    to :math:`\\phi`. A simple rule is to make it constant, for example
    :math:`\\phi = 0.5`. A much better approach due to Smock (1962), is to
    make :math:`\\phi = 1/n`.
    Parameters
    ----------
    graph : NetworkX graph
       Graph which represents the road network with the actual locations of
       the source and sink nodes and the connections to the road network.
    od_graph : NetworkX graph
       Graph contains a network between source and sink notes. The demand of
       each edge should be already assigned.
    od_matrix : numpy matrix
       Matrix with origin destination demands. If the od_graph does not
       contain this information, it can be optional loded as separate matrix.
    Examples
    --------
    >>> traffic = TrafficModel(graph,od_graph) # doctest: +SKIP
    >>> traffic.run()
    >>> graph_results = traffic.get_graph()
    References
    ----------
    .. [1] Juan de Dios Ortuzar and Luis G. Willumsen (2011) Modelling Transport, Fourth Edition. John Wiley and Sons.
    """

    def __init__(self, graph, od_graph, od_matrix=None, fcoeffs=[1,0,0,0,0.15,0], threshold=5e-4, iterations=400):

        self.graph = graph
        self.od_graph = od_graph
        self.od_matrix = od_matrix
        self.lost_trips = {}
        self.cut_links = []

        self.unit_factor = 1000
        self.fcoeffs = fcoeffs
        self.threshold = threshold
        self.large = 1e+14
        # number of iterations for the corrective factors
        self.n_iter_tm = iterations

        self.temp_list = []
        self.od_paths = collections.defaultdict(dict)
        pass

    def get_edge_attribute(self, edge, attribute):
        return self.graph[edge[0]][edge[1]][attribute]

    def set_edge_attribute(self, edge, attribute, value):
        self.graph[edge[0]][edge[1]][attribute] = value

    def calculate_initial_traveltime(self):
        [self.set_edge_attribute(edge, 't_0', self.get_edge_attribute(edge, 'length') \
                                 / self.unit_factor / self.get_edge_attribute(edge, 'speedlimit'))\
                                for edge in self.graph.edges()]

    def set_initial_traveltimes(self):
        [self.set_edge_attribute(edge, 't_k', self.get_edge_attribute(edge, 't_0')) for edge in self.graph.edges()]
        [self.set_edge_attribute(edge, 't_h', 0) for edge in self.graph.edges()]

    def set_initial_flow(self):
        [self.set_edge_attribute(edge, 'flow', 0) for edge in self.graph.edges()]

    def set_initial_help_flow(self):
        [self.set_edge_attribute(edge, 'help_flow', 0) for edge in self.graph.edges()]

    def set_od_matix(self, od_matrix):
        for s,t in self.od_graph.edges():
            self.od_graph[s][t]['demand'] = od_matrix[s, t]

    #@timeit
    def set_help_traveltime(self):
        [self.set_edge_attribute(edge, 't_h', \
                self.get_edge_attribute(edge, 't_k')) for edge in self.graph.edges()]


    def solve_and_assign_sp(self, ods, weight='t_k'):
        s,t = ods
        dist, sp = nx.single_source_dijkstra(self.graph, source=s, target=t, weight=weight)
        LBDst = dist * self.od_graph[s][t]['demand']
        dic = {}
        [dic.setdefault((sp[i],sp[i + 1]), self.od_graph[s][t]['demand']) for i in range(len(sp) - 1)]
        #for i in range(len(sp) - 1):
        #    dic[(sp[i],sp[i + 1])]=self.od_graph[s][t]['demand']
        return LBDst, dic

    #@timeit
    def calculate_auxiliary_flows(self, pool):
        self.set_initial_help_flow()
        self.LBD = 0
        ods = [(s,t) for s,t in self.od_graph.edges()]
        results = pool.map(self.solve_and_assign_sp, ods)
        self.LBD += sum([lbds for lbds, flows in results])
        for lbds, flows in results:
            for a, flow in flows.items():
                self.graph[a[0]][a[1]]['help_flow'] += flow
    #@timeit
    def calculate_flows(self, phi):
        [self.set_edge_attribute(edge, 'flow', (1 - phi) * self.get_edge_attribute(edge, 'flow') + \
            phi * self.get_edge_attribute(edge, 'help_flow')) or edge in self.graph.edges() for edge in self.graph.edges()]


    def calculate_traveltime_social_edge_exo(self, edge,  fcoeffs, exogenous_G=False):
         self.set_edge_attribute(edge, 't_k',  self.get_edge_attribute(edge, 't_0') * sum([fcoeffs[i] * \
                                    ((self.get_edge_attribute(edge, 'flow')+exogenous_G[edge[0]][edge[1]]['flow'])/self.get_edge_attribute(edge, 'capacity'))**(i) for i in range(len(fcoeffs))]) \
                                    + self.get_edge_attribute(edge, 't_0') *  sum([fcoeffs[i+1]*(i+1)\
                                    *((self.get_edge_attribute(edge, 'flow')+exogenous_G[edge[0]][edge[1]]['flow'])/self.get_edge_attribute(edge, 'capacity'))**(i) for i in range(len(fcoeffs)-1)]))

    def calculate_traveltime_social_edge(self, edge,  fcoeffs):
        self.set_edge_attribute(edge, 't_k',  self.get_edge_attribute(edge, 't_0') * sum([fcoeffs[i] * \
                                ((self.get_edge_attribute(edge, 'flow'))/self.get_edge_attribute(edge, 'capacity'))**(i) for i in range(len(fcoeffs))]) \
                                + self.get_edge_attribute(edge, 't_0') *  sum([fcoeffs[i+1]*(i+1)\
                                *((self.get_edge_attribute(edge, 'flow'))/self.get_edge_attribute(edge, 'capacity'))**(i) for i in range(len(fcoeffs)-1)]))

    #@timeit
    def calculate_traveltime_social(self, fcoeffs, exogenous_G=False):
        # TODO: add description
        if exogenous_G == False:
            [self.calculate_traveltime_social_edge(edge,  fcoeffs) for edge in self.graph.edges()]
        else:
            [self.calculate_traveltime_social_edge_exo(edge,  fcoeffs, exogenous_G=exogenous_G) for edge in self.graph.edges()]

    #@timeit
    def calculate_traveltime(self, fcoeffs, exogenous_G=False):
        if exogenous_G == False:
            [self.set_edge_attribute(edge, 't_k', self.get_edge_attribute(edge, 't_0') * \
                    sum([fcoeffs[i] * ((self.get_edge_attribute(edge, 'flow'))/self.get_edge_attribute(edge, 'capacity'))**(i) for i in range(len(fcoeffs))])) \
                    for edge in self.graph.edges()]
        else:
            [self.set_edge_attribute(edge, 't_k', self.get_edge_attribute(edge, 't_0') * \
                sum([fcoeffs[i] * ((self.get_edge_attribute(edge, 'flow'))+exogenous_G[edge[0]][edge[1]]['flow'] \
                                   / self.get_edge_attribute(edge,'capacity')) ** (i) for i in range(len(fcoeffs))])) \
                for edge in self.graph.edges()]

    #@timeit
    def stopping_criteria(self):
        TT = self.calculate_total_travel_time()
        self.RG = abs((TT-self.LBD)/self.LBD)
        if self.RG <= self.threshold:
            return True
        else:
            return  False

    def check_network_connections(self):
        graph = self.graph.copy()
        [graph.remove_edge(edge[0], edge[1]) for edge in graph.edges() if self.get_edge_attribute(edge, 'capacity') == 0]
        cut_links = []
        for s,t in list(self.od_graph.edges()):
            if not nx.has_path(graph, s, t):
                self.lost_trips[(s, t)] = self.od_graph[s][t]['demand']
                self.od_graph.remove_edge(s, t)
                cut_value, partition = nx.minimum_cut(self.graph, s, t)
                reachable, non_reachable = partition
                cutset = set()
                [cutset.update((u, v) for v in nbrs if v in non_reachable) for u, nbrs in ((n, self.graph[n]) for n in reachable)]
                cut_links.extend(list(cutset))
        [self.cut_links.append(edge) for edge in list(set(cut_links)) if self.graph[edge[0]][edge[1]]['capacity'] == 0]

    def calculate_total_travel_time(self):
        return sum([self.get_edge_attribute(edge, 't_k')*self.get_edge_attribute(edge, 'flow') for edge in self.graph.edges()])

    #@timeit
    def calculate_phi(self, bisection_n=10, exogenous_G=False):
        l = 1/2
        ub = 1
        lb = 0
        for n in range(bisection_n): #Solve issue with exogenous
            dx1dl = sum(get_travel_time(x=l*self.get_edge_attribute(edge, 'help_flow')+(1-l)*self.get_edge_attribute(edge, 'flow'),\
                                    m=self.get_edge_attribute(edge, 'capacity'),\
                                    t0=self.get_edge_attribute(edge, 't_0'),
                                    fcoeffs=self.fcoeffs) \
                    * (self.get_edge_attribute(edge, 'help_flow') - self.get_edge_attribute(edge, 'flow')) \
                    for edge in self.graph.edges())
            if dx1dl>0:
                ub = l
            else:
                lb = l
            l = (ub - lb) / 2
        return l

    #@timeit
    def calculate_phi_social(self, bisection_n=10, exogenous_G=False):
        l = 1/2
        ub = 1
        lb = 0
        for n in range(bisection_n): #Solve issue with exogenous
            dx1dl = sum(get_travel_time_social(x=l*self.get_edge_attribute(edge, 'help_flow')+(1-l)*self.get_edge_attribute(edge, 'flow'),\
                                    m=self.get_edge_attribute(edge, 'capacity'),\
                                    t0=self.get_edge_attribute(edge, 't_0'),
                                    fcoeffs=self.fcoeffs) \
                    * (self.get_edge_attribute(edge, 'help_flow') - self.get_edge_attribute(edge, 'flow')) \
                    for edge in self.graph.edges())
            if dx1dl>0:
                ub = l
            else:
                lb = l
            l = (ub - lb) / 2

        return l


    def run(self, fcoeffs=[1,0,0,0,0.15,0], build_t0=False, exogenous_G=False):
        pool = mp.Pool(mp.cpu_count()-1)
        # assign od matrix to od graph (if matrix is given)
        if self.od_matrix is not None:
            self.set_od_matix(self.od_matrix)
        # calculate traveltime at t=0
        if build_t0:
            self.calculate_initial_traveltime()
        # set traveltime equal to initial traveltime
        self.set_initial_traveltimes()
        self.assert_t0_exists()
        # set initial flow = 0
        self.set_initial_flow()
        # set initial help flow = 0
        self.set_initial_help_flow()
        # check network if every source and target can be reached
        self.check_network_connections()
        phi = 1
        self.RG=99
        for i in range(1, self.n_iter_tm):
            # calculate auxiliary flows
            self.calculate_auxiliary_flows(pool)
            # calculating the flow using auxiliary flow
            self.calculate_flows(phi)
            # save old traveltimes
            self.set_help_traveltime()
            # calculate new traveltime
            self.calculate_traveltime(fcoeffs=fcoeffs, exogenous_G=exogenous_G)
            #update phi
            phi = self.calculate_phi()
            # stopping criteria
            print("i: " + str(i) + '\t RG: ' + str(self.RG))
            if self.stopping_criteria():
                #print("i: " + str(i) + '\t RG: ' + str(self.RG))
                pool.close()
                break
        pool.close()

    def run_social(self, fcoeffs=[1, 0, 0, 0, 0.15, 0], build_t0=False, exogenous_G=False):
        pool = mp.Pool(mp.cpu_count()-1)
        # assign od matrix to od graph (if matrix is given)
        if self.od_matrix is not None:
            self.set_od_matix(self.od_matrix)
        # calculate traveltime at t=0
        if build_t0:
            self.calculate_initial_traveltime()
        # set traveltime equal to initial traveltime
        self.set_initial_traveltimes()
        self.assert_t0_exists()
        # set initial flow = 0
        self.set_initial_flow()
        # set initial help flow = 0
        self.set_initial_help_flow()
        # check network if every source and target can be reached
        self.check_network_connections()
        phi = 1
        self.RG = 99
        for i in range(1, self.n_iter_tm):
            # calculate auxiliary flows
            self.calculate_auxiliary_flows(pool)
            # calculating the flow using auxiliary flow
            self.calculate_flows(phi)
            # save old traveltimes
            self.set_help_traveltime()
            # calculate new traveltime
            self.calculate_traveltime_social(fcoeffs=fcoeffs, exogenous_G=exogenous_G)
            # update phi
            phi = self.calculate_phi()
            # stopping criteria
            print("i: " + str(i-1) + '\t RG: ' + str(self.RG))
            if self.stopping_criteria():
                print("i: " + str(i) + '\t RG: ' + str(self.RG))
                pool.close()
                break
        pool.close()



    def get_traveltime(self):
        traveltime = []
        for edge in self.graph.edges():
            if self.get_edge_attribute(edge, 'type') != 'Misc':
                if self.get_edge_attribute(edge, 'capacity') != 0:
                    traveltime.append(self.get_edge_attribute(edge, 't_k'))
        return traveltime

    def get_flow(self):
        flow = []
        for edge in self.graph.edges():
            if self.get_edge_attribute(edge, 'type') != 'Misc':
                if self.get_edge_attribute(edge, 'capacity') != 0:
                    flow.append(self.get_edge_attribute(edge, 'flow'))
        return flow

    def get_total_flow_estimated(self):
        return sum([sum(self.od_paths[od].values()) for od in self.od_paths])

    def get_flow_error(self):
        est = sum([sum(self.od_paths[od].values()) for od in self.od_paths])
        ini = np.sum(self.od_matrix)
        return est-ini, est/ini - 1

    def get_od_path_flows(self):
        return self.od_paths

    def get_car_hours(self):
        car_hours = []
        for edge in self.graph.edges():
            if self.get_edge_attribute(edge, 'type') != 'Misc':
                if self.get_edge_attribute(edge, 'capacity') != 0:
                    car_hour = self.get_edge_attribute(
                        edge, 't_k') * self.get_edge_attribute(edge, 'flow')
                    car_hours.append(car_hour)
        return car_hours

    def get_car_distances(self):
        """vehicle-kilometre"""
        car_distances = []
        for edge in self.graph.edges():
            if self.get_edge_attribute(edge, 'type') != 'Misc':
                if self.get_edge_attribute(edge, 'capacity') != 0:
                    car_distance = self.get_edge_attribute(
                        edge, 'length') / self.unit_factor * self.get_edge_attribute(edge, 'flow')
                    car_distances.append(car_distance)
        return car_distances

    def get_lost_trips(self):
        return self.lost_trips

    def get_cut_links(self):
        return self.cut_links

    def print_results(self):
        for edge in self.graph.edges():
            if self.get_edge_attribute(edge, 'type') != 'Misc':
                name = self.get_edge_attribute(edge, 'name')
                flow = self.get_edge_attribute(edge, 'flow')
                initialtraveltime = self.get_edge_attribute(edge, 't_0')
                traveltime = self.get_edge_attribute(edge, 't_k')
                print('t_0 vs t_k (flow) ' + str(name) + ': ' + str(round(initialtraveltime, 2)
                                                                    ) + ', ' + str(round(traveltime, 2)) + ' ('+str(round(flow, 0))+')')

    def get_graphs(self):
        return self.graph, self.od_graph

    def get_od_matrix(self):
        return self.od_matrix

    def get_graph(self):
        return self.graph

    def assert_t0_exists(self):
        for edge in self.graph.edges():
            assert (self.get_edge_attribute(edge,'t_0') > 0) ,"Edge: "+str(edge)+" has t_0 = 0 or t_0 not exists!!"



def get_travel_time(x, m, t0,  fcoeffs, exo=0):
        return t0 *sum(fcoeffs[i]*(x+exo/m)**i for i in range(len(fcoeffs)))

def get_travel_time_social(x, m, t0, fcoeffs, exo=0):
    return t0 *sum(fcoeffs[i]*(x+exo/m)**i for i in range(len(fcoeffs))) + t0 * sum(fcoeffs[i+1]*(i+1)*(x+exo/m)**i for i in range(len(fcoeffs)-1))


def get_dxdb(TAP, delta=0.05, divide=1, num_cores = False):
    """
    Derivative of the flow distribution with respect to all beta coefficients
    in the cost function. Uses a finite-difference approach

    parameters
    ----------

    fcoeffs: the actual coefficient vector
    delta: how big/small is the step

    Returns
    -------
    a list (position is the number of parameter of the polynomial) 
    with dictionaries (key=edge, value=derivative)
    """
    nPoly = len(TAP.fcoeffs)
    if num_cores == False:
        num_cores = multiprocessing.cpu_count()
    G_original = TAP.graph.copy()
    fcoeffs_i = []
    fcoeffs_i_j = TAP.fcoeffs.copy()
    for i in range(nPoly):  
        e = e_vect(nPoly, i)
        fcoeffs_i.append( [TAP.fcoeffs[coeff] + delta*e[coeff] for coeff in range(nPoly)] ) 
    #dxdb = Parallel(n_jobs=num_cores)(delayed(get_dxdb_sing)(TAP, i) for i in fcoeffs_i)
    pool = Pool(num_cores)
    dxdb = [pool.apply(get_dxdb_sing, (TAP, i)) for i in fcoeffs_i]
    pool.close()
    return dxdb



def get_dxdb_sing(TAP, fcoeffs_i, delta=0.05, divide=1):
    """
    Derivative of the flow distribution with respect to a single beta coefficient
    in the cost function. Uses a finite-difference approach

    parameters
    ----------

    fcoeffs: the actual coefficient vector
    delta: how big/small is the step

    Returns
    -------
    a dictionary (key=edge, value=derivative)
    """
    G_original = TAP.graph.copy()
    TAP.run(fcoeffs=fcoeffs_i, build_t0=False)
    G_modified = TAP.graph.copy()
    flowDiff_ = flowDiff(G_modified, G_original)
    dxdbi = {k : v/delta/divide for k,v in flowDiff_.items()}
    return dxdbi 

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

@timeit
def get_dxdg(G, gDict, k=1):
    """
    Derivative of the flow distribution with respect to the OD demands.

    Parameters
    ----------

    G: a network object with "flow" attribute
    gDict: It's corresponding OD demand

    Returns
    -------
    a sparse matrix with the edges as keys and derivative as value
    """
    
    path = [(s, {t: kShortestPaths(G, s, t, k, weight="t_k") for a,t in gDict.keys()} ) for s,b in gDict.keys()] 
    od_route_dict = {}
    for s,t in gDict.keys():
        route = str(path[s-1][1][t]).replace("[", "").replace(", ", "->").replace("]", "")
        od_route_dict[(s,t)] = route

    od_link_dict = {}
    for s,t in od_route_dict.keys():
        od_link_list = []
        od_node_list = od_route_dict[(s,t)].split('->')
        for i in range(len(od_node_list)):
            if i < len(od_node_list) - 1:
                od_link_list.append((int(od_node_list[i]), int(od_node_list[i+1])))
        od_link_dict[(s,t)] = od_link_list
    jacob = {}
    for s,t in gDict.keys():
        for u,v in G.edges():
            if (u,v) in od_link_dict[(s,t)]:
                jacob[(s,t), (u,v)] = 1
            else:
                jacob[(s,t), (u,v)] = 0
    return jacob




def kShortestPaths(G, source, target, k, weight=None):
    """
    Returns the k shortest paths from source to target

    Parameters
    ----------
    G: a networkx graph
    source: the starting node
    target: the goal node
    k: an integer of number of routes
    weight: if a weight is going to be considered

    Returns
    -------
    a list containing shortest path objects. 

    """
    return list(itertools.islice(nx.shortest_simple_paths(G, source, target, weight=weight), k))


# =============================================================================
# eof
#
# Local Variables:
# mode: python
# mode: linum
# End: