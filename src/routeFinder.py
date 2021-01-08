from  gurobipy import *
import networkx as nx
import numpy as np
from itertools import islice
import src.tnet as tnet
from src.utils import *

def k_shortest_paths(G, source, target, k, weight=None):
    return list(islice(nx.shortest_simple_paths(G, source, target, weight=weight), k))


def solve_flow_finder(G, gw, R, xw):
    # Start model
    m = Model('QP')
    m.setParam('OutputFlag',0)
    m.setParam('BarHomogeneous', 1)
    m.update()

    links = list(G.edges())
    nLinks = len(links)
    nRoutes = len(R)
    gw = gw[1]

    xw = [x for k,x in xw.items()]

    # Define variables
    p = [m.addVar(lb=0, ub=1, name='p'+str(i)) for i in range(nRoutes)]
    #x = [m.addVar(lb=0, ub=1, name='x_' + str(i) + '_'+str(j)) for (i,j) in G.edges()]
    m.update()

    # Build route-link incidence Matrix
    A = [[0 for i in range(nRoutes)] for j in range(nLinks)]
    ri = 0
    for r in R:
        for i in range(len(r)-1):
            l = links.index((r[i], r[i+1]))
            A[l][ri] = 1
        ri += 1

    A = np.array(A)

    # Set Objective
    x = []
    for l in range(nLinks):
        x.append(sum([A[l][r] * p[r] * gw for r in range(nRoutes)]))

    obj = sum([(xw[i]-x[i])*(xw[i]-x[i]) for i in range(len(x))])
    m.update()
    m.setObjective(obj, GRB.MINIMIZE)

    # Set Constaints
    m.addConstr(sum(p) == 1)
    m.update()
    m.optimize()
    sol = get_sol(m, R)
    return sol, obj.getValue()


def get_sol(m, R):
    sol_dic = {}
    j = 0
    for r in range(len(R)):
        p = m.getVarByName('p' + str(r)).X
        if p>0.001:
            sol_dic[j] = {}
            sol_dic[j]['p'] = p
            sol_dic[j]['r'] = R[r]
            j+=1
    return  sol_dic

@timeit
def routeFinder_OD(G, gw, od_flows_w, eps, max_routes=100):
    i, j = gw[0]
    k = 0
    err = 999999
    while (err >= eps) and (k<max_routes):
        k += 1
        R = k_shortest_paths(G, source=i, target=j, k=k, weight='t_k')
        sol, obj = solve_flow_finder(G, gw, R, od_flows_w)
        err = obj
        #print(str(k) + ' : ' + str(err))
    return sol


def routeFinder(G, g, od_flows,  eps):
    routes = {}
    for gi in g.items():
        routes[gi[0]] = routeFinder_OD(G, gi, od_flows[gi[0]], eps=eps)
    return routes


## CODE FOR REBALANCING FINDER


def get_node_potentials(G):
    for n in G.nodes():
        in_flow = sum([G[i][j]['flowRebalancing'] for i, j in G.in_edges(n)])
        out_flow = sum([G[i][j]['flowRebalancing'] for i, j in G.out_edges(n)])
        G[n]['potential'] = in_flow - out_flow

def get_rebalancing_ods(G, eps = 0.0):
    o = []
    d = []
    for n in G.nodes():
        in_flow = sum([G[i][j]['flowRebalancing'] for i,j in G.in_edges(n)])
        out_flow  =  sum([G[j][k]['flowRebalancing'] for j,k in G.out_edges(n)])
        G.nodes[n]['potential'] = in_flow - out_flow
        if in_flow - out_flow >= eps:
            d.append(n)
        elif in_flow - out_flow <= -eps:
            o.append(n)
    return o,d

def get_sorted_ods(G, o, d, astar=False):
    r = {origin:{dest: nx.shortest_path_length(G, origin, dest, 't_0') for dest in d} for origin in o}
    o_list = {}
    for origin in o:
        o_list[origin] = [k for k,v in sorted(r[origin].items(), key=lambda item: item[1])]
    d_list = {}
    for dest in d:
        d_list[dest] = [a[1] for a in sorted([[r[origin][dest],origin] for origin in o])]
    return o_list, d_list

@timeit
def solve_rebalancing_O(G, O):
    m = Model()
    m.setParam('OutputFlag', 0)
    m.setParam('BarHomogeneous', 1)
    m.setParam('Method', 1)
    # Define vars
    [m.addVar(lb=0, name='x^'+str(o)+'_'+str(i)+'_'+str(j)) for i,j in G.edges() for o in O]
    m.update()
    x = {}
    for i,j in G.edges():
        x[(i,j)] = quicksum([m.getVarByName('x^'+str(o)+'_'+str(i)+'_'+str(j)) for o in O])
    # Add Obj
    obj = quicksum(G[i][j]['t_k']*x[(i,j)] for i,j in G.edges())
    m.setObjective(obj)
    m.update()
    # Add constraints
    [m.addConstr(quicksum(x[(i,j)] for i,j in G.in_edges(nbunch=n)) - quicksum(x[(j,k)] for j,k in G.out_edges(nbunch=n)) == G.nodes[n]['potential']) for n in G.nodes()]
    [m.addConstr(x[(i, j)]-G[i][j]['flowRebalancing'] == 0) for i,j in G.edges()]
    m.update()
    for o in O:
        for n in G.nodes() :
            if n == o :
                m.addConstr(quicksum(m.getVarByName('x^' + str(o) + '_' + str(i) + '_' + str(j)) for i, j in G.in_edges(nbunch=n))
                            - quicksum(m.getVarByName('x^' + str(o) + '_' + str(j) + '_' + str(k)) for j, k in G.out_edges(nbunch=n)) == G.nodes[n]['potential'])
            else:
                m.addConstr(quicksum(m.getVarByName('x^' + str(o) + '_' + str(i) + '_' + str(j)) for i, j in G.in_edges(nbunch=n))
                            - quicksum(m.getVarByName('x^' + str(o) + '_' + str(j) + '_' + str(k)) for j, k in G.out_edges(nbunch=n)) >= 0)
    m.update()
    m.optimize()
    return {o:{(i,j): m.getVarByName('x^'+str(o)+'_'+str(i)+'_'+str(j)).X for i,j in G.edges()} for o in O}

@timeit
def solve_rebalancing_D(G, origin, D, xo, potential):
    m = Model()
    m.setParam('OutputFlag', 0)
    m.setParam('Method', 1)
    # Define vars
    [m.addVar(lb=0, name='x^'+str(d)+'_'+str(i)+'_'+str(j)) for i,j in G.edges() for d in D]
    m.update()
    x = {}
    for i,j in G.edges():
        x[(i,j)] = quicksum([m.getVarByName('x^'+str(d)+'_'+str(i)+'_'+str(j)) for d in D])
    m.update()
    # Add Obj
    obj = quicksum(G[i][j]['t_k']*x[(i,j)] for i,j in G.edges())
    m.setObjective(obj)
    m.update()
    # Add constraints
    [m.addConstr(x[(i, j)] - xo[(i, j)] == 0) for i, j in G.edges()]
    m.addConstr(quicksum(x[(i,j)] for i,j in G.in_edges(nbunch=origin)) - quicksum(x[(j,k)] for j,k in G.out_edges(nbunch=origin)) == potential[origin])
    m.update()
    for d in D:
        for n in G.nodes() :
            if n == origin:
                m.addConstr(quicksum(m.getVarByName('x^' + str(d) + '_' + str(i) + '_' + str(j)) for i, j in G.in_edges(nbunch=n))
                            - quicksum(m.getVarByName('x^' + str(d) + '_' + str(j) + '_' + str(k)) for j, k in G.out_edges(nbunch=n)) <= 0)
            elif n == d:
                m.addConstr(quicksum(m.getVarByName('x^' + str(d) + '_' + str(i) + '_' + str(j)) for i, j in G.in_edges(nbunch=n))
                            - quicksum(m.getVarByName('x^' + str(d) + '_' + str(j) + '_' + str(k)) for j, k in G.out_edges(nbunch=n)) >= 0)
            else:
                m.addConstr(quicksum(m.getVarByName('x^' + str(d) + '_' + str(i) + '_' + str(j)) for i, j in G.in_edges(nbunch=n))
                            - quicksum(m.getVarByName('x^' + str(d) + '_' + str(j) + '_' + str(k)) for j, k in G.out_edges(nbunch=n)) == 0)

    m.update()
    m.optimize()
    m.update()
    return {d: {(i, j): m.getVarByName('x^' + str(d) + '_' + str(i) + '_' + str(j)).X for i, j in G.edges()} for d in D}



@timeit
def solve_flow_decomposition_D(G, origin, xo, g, L=1):
    m = Model()
    m.setParam('OutputFlag', 0)
    m.setParam('Method', 1)
    # Define vars
    D = [d for o,d in g.keys() if o==origin]
    [m.addVar(lb=0, name='x^'+str(d)+'_'+str(i)+'_'+str(j)) for i,j in G.edges() for d in D]
    m.update()
    x = {}
    for i,j in G.edges():
        x[(i,j)] = quicksum([m.getVarByName('x^'+str(d)+'_'+str(i)+'_'+str(j)) for d in D])
    m.update()
    # Add Obj
    obj = quicksum(L*G[i][j]['t_k']*x[(i,j)] for i,j in G.edges())
    m.setObjective(obj)
    m.update()
    # Add constraints
    [m.addConstr(x[(i, j)] - xo[(i, j)] == 0) for i, j in G.edges()]
    #m.addConstr(quicksum(x[(i,j)] for i,j in G.in_edges(nbunch=origin)) - quicksum(x[(j,k)] for j,k in G.out_edges(nbunch=origin)) == potential[origin])
    m.update()
    for d in D:
        for n in G.nodes() :
            if origin != d:
                if n == origin:
                    m.addConstr(quicksum(m.getVarByName('x^' + str(d) + '_' + str(i) + '_' + str(j)) for i, j in G.in_edges(nbunch=n))
                                 - quicksum(m.getVarByName('x^' + str(d) + '_' + str(j) + '_' + str(k)) for j, k in G.out_edges(nbunch=n)) + g[(origin, d)] == 0)
                                #- quicksum(m.getVarByName('x^' + str(d) + '_' + str(j) + '_' + str(k)) for j, k in G.out_edges(nbunch=n)) <= 0)
                elif n == d:
                    m.addConstr(quicksum(m.getVarByName('x^' + str(d) + '_' + str(i) + '_' + str(j)) for i, j in G.in_edges(nbunch=n))
                                - g[(origin, d)]- quicksum(m.getVarByName('x^' + str(d) + '_' + str(j) + '_' + str(k)) for j, k in G.out_edges(nbunch=n)) ==0)
                else:
                    m.addConstr(quicksum(m.getVarByName('x^' + str(d) + '_' + str(i) + '_' + str(j)) for i, j in G.in_edges(nbunch=n))
                                - quicksum(m.getVarByName('x^' + str(d) + '_' + str(j) + '_' + str(k)) for j, k in G.out_edges(nbunch=n))==0)
            else:
                m.addConstr(quicksum(m.getVarByName('x^' + str(d) + '_' + str(i) + '_' + str(j)) for i, j in G.in_edges(nbunch=n))
                            - quicksum(m.getVarByName('x^' + str(d) + '_' + str(j) + '_' + str(k)) for j, k in G.out_edges(nbunch=n)) == 0)

    m.update()
    m.optimize()
    m.update()
    return {d: {(i, j): m.getVarByName('x^' + str(d) + '_' + str(i) + '_' + str(j)).X for i, j in G.edges()} for d in D}



@timeit
def solve_decomposition_dest(G, D, xo, potential, origin, eps=0.001):
    m = Model()
    m.setParam('OutputFlag', 0)
    m.setParam('Method', 1)
    # Define vars
    [m.addVar(lb=0, name='x^'+str(d)+'_'+str(i)+'_'+str(j)) for i,j in G.edges() for d in D]
    m.update()
    x = {}
    for i,j in G.edges():
        x[(i,j)] = quicksum([m.getVarByName('x^'+str(d)+'_'+str(i)+'_'+str(j)) for d in D])
    m.update()
    # Add Obj
    obj = quicksum(G[i][j]['t_k']*x[(i,j)] for i,j in G.edges())
    m.setObjective(obj)
    m.update()
    # Add constraints
    [m.addConstr(x[(i, j)] - xo[(i, j)] >= -eps) for i, j in G.edges()]
    [m.addConstr(x[(i, j)] - xo[(i, j)] <= eps) for i, j in G.edges()]

    for d in D:
        for j in G.nodes():
            if j == origin:
                m.addConstr(quicksum(m.getVarByName('x^' + str(d) + '_' + str(i) + '_' + str(origin)) for i, j in G.in_edges(nbunch=origin))
                - quicksum(m.getVarByName('x^' + str(d) + '_' + str(origin) + '_' + str(k)) for j, k in G.out_edges(nbunch=origin))
                == -potential[d])
            elif j ==d:
                m.addConstr(quicksum(m.getVarByName('x^' + str(d) + '_' + str(i) + '_' + str(j)) for i, j in G.in_edges(nbunch=j))
                            - quicksum(m.getVarByName('x^' + str(d) + '_' + str(j) + '_' + str(k)) for j, k in G.out_edges(nbunch=j))
                            == potential[d])
            else:
                m.addConstr(quicksum(m.getVarByName('x^' + str(d) + '_' + str(i) + '_' + str(j)) for i, j in G.in_edges(nbunch=j))
                            - quicksum(m.getVarByName('x^' + str(d) + '_' + str(j) + '_' + str(k)) for j, k in G.out_edges(nbunch=j))
                            <= 0)
    m.update()
    m.optimize()
    m.update()
    return {d: {(i, j): m.getVarByName('x^' + str(d) + '_' + str(i) + '_' + str(j)).X for i, j in G.edges()} for d in D}



@timeit
def get_destinations(G,x, eps=0.01):
    potential = {}
    o = []
    d = []
    for n in G.nodes():
        in_flow = sum([x[(i,j)] for i, j in G.in_edges(n)])
        out_flow = sum([x[(j,k)] for j, k in G.out_edges(n)])
        potent = in_flow - out_flow
        if potent >= eps:
            potential[n] = potent
            d.append(n)
        elif potent <= -eps:
            o.append(n)
            potential[n] = potent
        else:
            potential[n] = 0
    return o, d, potential

@timeit
def rebRouteFinder(G, eps, print_=False):
    routes_dic = {}
    O,D = get_rebalancing_ods(G, eps=0)
    xo = solve_rebalancing_O(G, O)
    if print_:
        print('Origin-flows have been found!')
    for origin, x in xo.items():
        Oo, Do, potential = get_destinations(G,x, eps=0.001)
        #xf = solve_rebalancing_D(G, origin, Do, x, potential)
        xf = solve_decomposition_dest(G, Do, x, potential, origin)
        if print_:
            print('OD-flows have been found!')
        xf = {(origin, d):v for d, v in xf.items()}
        for o,d in xf.keys():
            gw = [(o,d), sum([v for k,v in xf[(o,d)].items() if k[0]==o])]
            if gw[1] > eps:
                routes = routeFinder_OD(G, gw, xf[(o,d)], eps=eps, max_routes=20)
                routes_dic[(o,d)] = routes
    if print_:
        print('Route-flows have been found!')
    return routes_dic



def userRouteFinder(G, g, s_flows, eps):
    routes_dic = {}
    for origin, x in s_flows.items():
        xf = solve_flow_decomposition_D(G, origin, x, g, L=1)
        xf = {(origin, d):v for d, v in xf.items()}
        for o,d in xf.keys():
            gw = [(o,d), sum([v for k,v in xf[(o,d)].items() if k[0]==o])]
            if gw[1] > eps:
                routes = routeFinder_OD(G, gw, xf[(o,d)], eps=20, max_routes=20)
                routes_dic[(o,d)] = routes
    return routes_dic


#xf = solve_flow_decomposition_D(G, origin, x, g, L=1)

def RouteFinder(G, g, s_flows, eps, od):
    o,d = od
    x = s_flows[o]
    xf = solve_flow_decomposition_D(G, o, x, g, L=1)
    xf = {(o, dest):v for dest, v in xf.items()}
    gw = [(o,d), sum([v for k,v in xf[(o,d)].items() if k[0]==o])]
    routes = routeFinder_OD(G, gw, xf[(o,d)], eps=eps, max_routes=10)
    return routes

def RouteTravelTime(G,path):
    return sum([G[path[i]][path[i+1]]['t_k'] for i in range(len(path)-1)])