
from gurobipy import *
import pwlf as pw
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from networkx.readwrite import json_graph
import json

def get_theta(fun):
    return fun.fit_breaks[1:-1]

def get_beta(fun):
    return fun.slopes[1:]

def get_approx_fun(fcoeffs, xa=False, range_=[0,2], nlines=3):
    x = list(np.linspace(range_[0], range_[1], 100))
    y = [sum( fcoeffs[n]*i**n for n in range(len(fcoeffs))) for i in x]

    my_pw = pw.PiecewiseLinFit(x, y)

    if xa == False:
        # fit the data for four line segments
        res = my_pw.fit(nlines)
    else:
        res = my_pw.fit_with_breaks([range_[0], xa,  range_[-1:][0]])




    # predict for the determined points
    xHat = np.linspace(min(x), max(x), num=10000)
    yHat = my_pw.predict(xHat)

    '''
    #plot the results
    plt.figure()
    plt.plot(x, y, '-')
    plt.plot(xHat, yHat, '-')
    plt.show()
    '''
    return  my_pw


def add_demand_cnstr(m, tnet):
    # Set Constraints
    for j in tnet.G_supergraph.nodes():
        for w, d in tnet.g.items():
            if j == w[0]:
                m.addConstr(quicksum(m.getVarByName('x^'+str(w)+'_'+str(i)+'_'+str(j)) for i,l in tnet.G_supergraph.in_edges(nbunch=j)) + d == quicksum(m.getVarByName('x^'+str(w)+'_'+str(j)+'_'+str(k)) for l,k in tnet.G_supergraph.out_edges(nbunch=j)))
            elif j == w[1]:
                m.addConstr(quicksum(m.getVarByName('x^' + str(w) + '_' + str(i) + '_' + str(j)) for i,l in tnet.G_supergraph.in_edges(nbunch=j)) == quicksum(m.getVarByName('x^' + str(w) + '_' + str(j) + '_' + str(k)) for l,k in tnet.G_supergraph.out_edges(nbunch=j)) + d)
            else:
                m.addConstr(quicksum(m.getVarByName('x^' + str(w) + '_' + str(i) + '_' + str(j)) for i,l in tnet.G_supergraph.in_edges(nbunch=j)) == quicksum(m.getVarByName('x^' + str(w) + '_' + str(j) + '_' + str(k)) for l,k in tnet.G_supergraph.out_edges(nbunch=j)))

    m.update()

def add_rebalancing_cnstr(m, tnet):
    for j in tnet.G.nodes():
        m.addConstr(quicksum(m.getVarByName('x^R'+str(i)+'_'+str(j)) + quicksum(m.getVarByName('x^' + str(w) + '_' + str(i) + '_' + str(j)) for w in tnet.g.keys() ) for i,l in tnet.G.in_edges(nbunch=j)) \
                    == quicksum(m.getVarByName('x^R'+str(j)+'_'+str(k)) + quicksum(m.getVarByName('x^' + str(w) + '_' + str(j) + '_' + str(k)) for w in tnet.g.keys() ) for l,k in tnet.G.out_edges(nbunch=j)) )
    m.update()

def set_optimal_flows(m,tnet):
    for i,j in tnet.G_supergraph.edges():
        if isinstance(i,str) or isinstance(j,str):
            tnet.G_supergraph[i][j]['flow'] = sum(m.getVarByName('x^' + str(w) + '_' + str(i) + '_' + str(j)).X for w in tnet.g.keys())
            tnet.G_supergraph[i][j]['flowNoRebalancing'] = sum(m.getVarByName('x^' + str(w) + '_' + str(i) + '_' + str(j)).X for w in tnet.g.keys())
        else:
            tnet.G_supergraph[i][j]['flow'] = m.getVarByName('x^R' + str(i) + '_' + str(j)).X + sum(m.getVarByName('x^' + str(w) + '_' + str(i) + '_' + str(j)).X for w in tnet.g.keys())
            tnet.G_supergraph[i][j]['flowNoRebalancing'] = sum(m.getVarByName('x^' + str(w) + '_' + str(i) + '_' + str(j)).X for w in tnet.g.keys())


def set_optimal_rebalancing_flows(m,tnet):
    for i,j in tnet.G.edges():
        tnet.G_supergraph[i][j]['flow'] += m.getVarByName('x^R' + str(i) + '_' + str(j)).X
        tnet.G_supergraph[i][j]['flowRebalancing'] = m.getVarByName('x^R' + str(i) + '_' + str(j)).X
        tnet.G_supergraph[i][j]['flowNoRebalancing'] = tnet.G[i][j]['flow'] - tnet.G[i][j]['flowRebalancing']


def set_CARS_par(tnet):
    # Set obj func parameters
    Vt = 100
    Vd = 0.1
    Ve = 1
    # Set the electricity constant
    ro = 1
    Af = 1
    cd = 1
    cr = 1
    mv = 1
    g = 9.81
    nu = 1
    for i,j in tnet.G.edges():
        tnet.G[i][j]['e'] =  (ro/2 *Af*cd * (tnet.G[i][j]['t_0']/tnet.G[i][j]['length'])**2 *cr * mv * g)* tnet.G[i][j]['length']/nu
    return Vt, Vd, Ve

def solve_CARS2(tnet, exogenous_G, fcoeffs, rebalancing=True, xa=0):
    #TODO: add description

    fun = get_approx_fun(fcoeffs, range_=[-.5,1.5])
    beta = get_beta(fun)
    theta = get_theta(fun)
    theta = [0.7, 1.2]
    beta = [0.53, 1.88]
    print(theta)
    print(beta)

    Vt, Vd, Ve = set_CARS_par(tnet)

    # Set exogenous flow
    if exogenous_G == 0:
        exogenous_G = tnet.G.copy()
        for i,j in tnet.G.edges():
            exogenous_G[i][j]['flow'] = 0

    # Start model
    m = Model('QP')
    m.setParam('OutputFlag',0)
    m.setParam('BarHomogeneous', 1)
    m.update()

    # Define variables
    var_xw = [m.addVar(lb=0, name='x^'+str(w)+'_'+str(i)+'_'+str(j)) for i,j in tnet.G_supergraph.edges() for w in tnet.g.keys()]
    var_xr = [m.addVar(lb=0, name='x^R'+str(i)+'_'+str(j)) for i,j in tnet.G.edges()]
    var_e1 = [m.addVar(lb=0, name='e^1_'+str(i)+'_'+str(j)) for i,j in tnet.G.edges()]
    var_e2 = [m.addVar(lb=0, name='e^2_'+str(i)+'_'+str(j)) for i,j in tnet.G.edges()]
    m.update()

    # Set objective
    obj = 0
    for i,j in tnet.G_supergraph.edges():
        for w in tnet.g.keys():
            obj += Vt *tnet.G_supergraph[i][j]['t_0'] * m.getVarByName('x^'+str(w)+'_'+str(i)+'_'+str(j))
    for i,j in tnet.G.edges():
        obj +=  Vt * (tnet.G[i][j]['t_0'] * beta[0]/tnet.G[i][j]['capacity']) * m.getVarByName('e^1_'+str(i)+'_'+str(j)) \
            * (m.getVarByName('e^1_'+str(i)+'_'+str(j)) + theta[0]*tnet.G[i][j]['capacity'] - exogenous_G[i][j]['flow']) \
            + Vt * (tnet.G[i][j]['t_0'] * beta[1]/tnet.G[i][j]['capacity']) * m.getVarByName('e^2_'+str(i)+'_'+str(j)) \
            * (m.getVarByName('e^2_'+str(i)+'_'+str(j)) + theta[1]*tnet.G[i][j]['capacity'] - exogenous_G[i][j]['flow']) \
            + Vt * (tnet.G[i][j]['t_0'] * beta[0] * m.getVarByName('e^2_'+str(i)+'_'+str(j)) * (theta[1]*tnet.G[i][j]['capacity'] - theta[0]*tnet.G[i][j]['capacity']) ) \
            + (Vd * tnet.G[i][j]['length'] + Ve * tnet.G[i][j]['e']) * (sum(
                m.getVarByName('x^' + str(w) + '_' + str(i) + '_' + str(j)) for w in tnet.g.keys()) + m.getVarByName(
                'x^R' + str(i) + '_' + str(j)))
    m.update()

    # Set Constraints
    add_demand_cnstr(m, tnet)
    if rebalancing==True:
        add_rebalancing_cnstr(m,tnet)
    for i,j in tnet.G.edges():
        m.addConstr(m.getVarByName('e^1_'+str(i)+'_'+str(j)) >= quicksum(m.getVarByName('x^' + str(w) + '_' + str(i) + '_' + str(j)) for w in tnet.g.keys()) +  m.getVarByName('x^R'+str(i)+'_'+str(j)) + exogenous_G[i][j]['flow'] - theta[0]*tnet.G[i][j]['capacity'] - m.getVarByName('e^2_'+str(i)+'_'+str(j)))
        m.addConstr(m.getVarByName('e^2_'+str(i)+'_'+str(j)) >= quicksum(m.getVarByName('x^' + str(w) + '_' + str(i) + '_' + str(j)) for w in tnet.g.keys()) +  m.getVarByName('x^R'+str(i)+'_'+str(j)) + exogenous_G[i][j]['flow'] - theta[1]*tnet.G[i][j]['capacity'])
    m.update()

    # Solve problem
    m.setObjective(obj, GRB.MINIMIZE)
    m.update()
    m.optimize()

    # saving  results
    set_optimal_flows(m, tnet)
    tnet.cars_obj = obj.getValue()
    return tnet


def solve_CARS(tnet, exogenous_G, fcoeffs, rebalancing=True, xa=0.01):
    fun = get_approx_fun(fcoeffs, xa=xa, nlines=2, range_=[0,1.5])
    beta = get_beta(fun)
    theta = get_theta(fun)

    #print(beta)
    #print(theta)

    Vt, Vd, Ve = set_CARS_par(tnet)
    # Set exogenous flow
    if exogenous_G == 0:
        exogenous_G = tnet.G.copy()
        for i,j in tnet.G.edges():
            exogenous_G[i][j]['flow'] = 0

    m = Model('QP')
    m.setParam('OutputFlag',0)
    m.setParam('BarHomogeneous', 1)
    m.update()

    # Define variables
    var_xw = [m.addVar(lb=0, name='x^'+str(w)+'_'+str(i)+'_'+str(j)) for i,j in tnet.G_supergraph.edges() for w in tnet.g.keys()]
    var_xr = [m.addVar(lb=0, name='x^R'+str(i)+'_'+str(j)) for i,j in tnet.G.edges()]
    var_e1 = [m.addVar(lb=0, name='e^1_'+str(i)+'_'+str(j)) for i,j in tnet.G.edges()]
    m.update()

    # Set objective
    obj = 0
    for i,j in tnet.G_supergraph.edges():
        for w in tnet.g.keys():
            obj += Vt * tnet.G_supergraph[i][j]['t_0'] * m.getVarByName('x^'+str(w)+'_'+str(i)+'_'+str(j))
    for i,j in tnet.G.edges():
        obj +=  Vt * (tnet.G[i][j]['t_0'] * beta[0]/tnet.G[i][j]['capacity']) * m.getVarByName('e^1_'+str(i)+'_'+str(j)) \
            * (m.getVarByName('e^1_'+str(i)+'_'+str(j)) + theta[0]*tnet.G[i][j]['capacity'] - exogenous_G[i][j]['flow']) \
            + (Vd * tnet.G[i][j]['length'] + Ve * tnet.G[i][j]['e']) * (sum(m.getVarByName('x^'+str(w)+'_'+str(i)+'_'+str(j)) for w in tnet.g.keys()) + m.getVarByName('x^R'+str(i)+'_'+str(j)))
    m.update()

    # Set Constraints
    add_demand_cnstr(m, tnet)
    if rebalancing==True:
        add_rebalancing_cnstr(m,tnet)
    for i,j in tnet.G.edges():
        m.addConstr(m.getVarByName('e^1_'+str(i)+'_'+str(j))  >= quicksum(m.getVarByName('x^' + str(w) + '_' + str(i) + '_' + str(j)) for w in tnet.g.keys()) +  m.getVarByName('x^R'+str(i)+'_'+str(j)) + exogenous_G[i][j]['flow'] - theta[0]* tnet.G[i][j]['capacity'])
    m.update()

    m.setObjective(obj, GRB.MINIMIZE)
    m.update()
    m.optimize()
    # saving  results
    set_optimal_flows(m, tnet)
    tnet.cars_obj = obj.getValue()
    return tnet


def solve_rebalancing(tnet, exogenous_G=0):
    Vt, Vd, Ve = set_CARS_par(tnet)
    # Set exogenous flow
    if exogenous_G == 0:
        exogenous_G = tnet.G.copy()
        for i, j in tnet.G.edges():
            exogenous_G[i][j]['flow'] = 0

    m = Model('QP')
    m.setParam('OutputFlag', 0)
    m.setParam('BarHomogeneous', 1)
    m.update()

    # Define variables
    var_xr = [m.addVar(lb=0, name='x^R' + str(i) + '_' + str(j)) for i, j in tnet.G.edges()]
    m.update()

    # Set objective
    obj = quicksum((Vd * tnet.G[i][j]['length'] + Ve * tnet.G[i][j]['e']) * m.getVarByName('x^R' + str(i) + '_' + str(j)) for i,j in tnet.G.edges())
    m.update()

    # Set Constraints
    for j in tnet.G.nodes():
        m.addConstr(quicksum(m.getVarByName('x^R' + str(i) + '_' + str(j)) + tnet.G[i][l]['flow'] for i, l in
                             tnet.G.in_edges(nbunch=j)) \
                    == quicksum(m.getVarByName('x^R' + str(j) + '_' + str(k)) + tnet.G[j][k]['flow'] for l, k in
                                tnet.G.out_edges(nbunch=j)))
    m.update()
    m.update()

    m.setObjective(obj, GRB.MINIMIZE)
    m.update()
    m.optimize()
    # saving  results
    set_optimal_rebalancing_flows(m,tnet)

def get_totalTravelTime_approx(tnet, fcoeffs, xa):
    fun = get_approx_fun(fcoeffs, xa=xa, nlines=2)
    beta = get_beta(fun)
    theta = get_theta(fun)
    print(beta)
    obj=0
    for i,j in tnet.G_supergraph.edges():
        if tnet.G_supergraph[i][j]['flow']/tnet.G_supergraph[i][j]['capacity'] <= xa:
            obj += tnet.G_supergraph[i][j]['flow'] * tnet.G_supergraph[i][j]['t_0']
        else:
            obj += tnet.G_supergraph[i][j]['flow'] * \
                   (tnet.G_supergraph[i][j]['t_0'] + (beta[0] *tnet.G_supergraph[i][j]['flow'] /tnet.G_supergraph[i][j]['capacity']))
    return obj


def travel_time(tnet, i, j, G_exo=False):
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
    if G_exo == False:
        return sum(
            [tnet.fcoeffs[n] * (tnet.G_supergraph[i][j]['flow'] / tnet.G_supergraph[i][j]['capacity']) ** n for n in
             range(len(tnet.fcoeffs))])
    else:
        return sum([tnet.fcoeffs[n] * ((tnet.G_supergraph[i][j]['flow'] + G_exo[i][j]['flow'])/ tnet.G_supergraph[i][j]['capacity']) ** n for n in range(len(tnet.fcoeffs))])


def get_totalTravelTime(tnet, G_exogenous=False):
    """
    evalute the travel time function on the SuperGraph level

    Parameters
    ----------

    tnet: transportation network object

    Returns
    -------
    float

    """
    if G_exogenous == False:
        return sum([tnet.G_supergraph[i][j]['flow'] * tnet.G_supergraph[i][j]['t_0'] * travel_time(tnet, i, j) for i, j in tnet.G_supergraph.edges()])
    else:
        ret = 0
        for i,j in tnet.G_supergraph.edges():
            if isinstance(tnet.G_supergraph[i][j]['type'], float)==True:
                ret += tnet.G_supergraph[i][j]['flow'] * tnet.G_supergraph[i][j]['t_0'] * travel_time_without_Rebalancing(tnet, i, j, G_exogenous[i][j]['flow'])
            else:
                ret += tnet.G_supergraph[i][j]['flow'] * tnet.G_supergraph[i][j]['t_0'] * travel_time_without_Rebalancing(tnet, i, j)
        return ret


def travel_time_without_Rebalancing(tnet, i, j, exo=0):
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
    return sum(
        [tnet.fcoeffs[n] * ((tnet.G_supergraph[i][j]['flow'] +exo )/ tnet.G_supergraph[i][j]['capacity']) ** n for n in range(len(tnet.fcoeffs))])

def get_totalTravelTime_without_Rebalancing(tnet, G_exogenous=False):
    """
    evalute the travel time function on the SuperGraph level

    Parameters
    ----------

    tnet: transportation network object

    Returns
    -------
    float

    """
    if G_exogenous==False:
        return sum([tnet.G_supergraph[i][j]['flowNoRebalancing'] * tnet.G_supergraph[i][j][
            't_0'] * travel_time_without_Rebalancing(tnet, i, j) for i, j in
                    tnet.G_supergraph.edges()])
    else:
        ret = 0
        for i,j in tnet.G_supergraph.edges():
            if isinstance(tnet.G_supergraph[i][j]['type'], float)==True:
                ret += tnet.G_supergraph[i][j]['flowNoRebalancing'] * tnet.G_supergraph[i][j]['t_0'] * travel_time_without_Rebalancing(tnet, i, j, G_exogenous[i][j]['flow'])
            else:
                ret += tnet.G_supergraph[i][j]['flowNoRebalancing'] * tnet.G_supergraph[i][j]['t_0'] * travel_time_without_Rebalancing(tnet, i, j)
    return ret



def get_pedestrian_flow(tnet):
    """
    get pedestrian flow in a supergraph

    Parameters
    ----------

    tnet: transportation network object

    Returns
    -------
    float

    """
    return sum([tnet.G_supergraph[i][j]['flow'] for i,j in tnet.G_supergraph.edges() if tnet.G_supergraph[i][j]['type']=='p'])

def get_amod_flow(tnet):
    """
    get amod flow in a supergraph

    Parameters
    ----------

    tnet: transportation network object

    Returns
    -------
    float

    """
    return sum([tnet.G_supergraph[i][j]['flowNoRebalancing'] for i,j in tnet.G.edges()])

def get_rebalancing_flow(tnet):
    """
    get rebalancing flow in a supergraph

    Parameters
    ----------

    tnet: transportation network object

    Returns
    -------
    float

    """
    return sum([tnet.G_supergraph[i][j]['flow']-tnet.G_supergraph[i][j]['flowNoRebalancing'] for i,j in tnet.G.edges()])

def plot_supergraph_car_flows(tnet, weight='flow', width=3, cmap=plt.cm.Blues):
    #TODO: add explaination
    fig, ax = plt.subplots()
    pos = nx.get_node_attributes(tnet.G, 'pos')
    d = {(i,j): tnet.G_supergraph[i][j][weight] for i,j in tnet.G.edges()}
    edges, weights = zip(*d.items())
    labels =  {(i,j): int(tnet.G_supergraph[i][j][weight]) for i,j in tnet.G.edges()}
    nx.draw(tnet.G, pos, node_color='b', edgelist=edges, edge_color=weights, width=width, edge_cmap=cmap)
    nx.draw_networkx_edge_labels(tnet.G, pos=pos, edge_labels=labels)
    return fig, ax

def plot_supergraph_pedestrian_flows(G, weight='flow', width=3, cmap=plt.cm.Blues):
	#TODO: add explaination
	fig, ax = plt.subplots()
	pos = nx.get_node_attributes(G, 'pos')
	edges, weights = zip(*nx.get_edge_attributes(G, weight).items())
	nx.draw(G, pos, node_color='b', edgelist=edges, edge_color=weights, width=width, edge_cmap=cmap)
	return fig, ax

def supergraph2G(tnet):
    # TODO: add explaination
    tnet.G = tnet.G_supergraph.subgraph(list([i for i in tnet.G.nodes()]))

def G2supergraph(tnet):
    # TODO: add explaination
    tnet.G_supergraph = tnet.G
    #for i,j in tnet.G_supergraph:
    #    tnet.G_supergraph[i][j]['flowNoRebalancing'] = tnet.G_supergraph[i][j]['flow']

def add_G_flows_no_rebalancing(array):
    # TODO: add description

    G = array[0].copy()
    for tn in array[1:]:
        for i, j in G.edges():
            G[i][j]['flow'] += tn[i][j]['flowNoRebalancing']
    return G


def solveMSAsocialCARS(tnet, exogenous_G=False):
    if exogenous_G == False:
        tnet.solveMSAsocial()
    else:
        tnet.solveMSAsocial(exogenous_G=exogenous_G)
    G2supergraph(tnet)


def nx2json(G, fname):
    with open(fname, 'w') as outfile1:
        outfile1.write(json.dumps(json_graph.node_link_data(G)))

def solve_social_Julia(tnet, exogenous_G=False):

    # Save to json files
    nx.write_gml(tnet.G, "tmp/G.txt")
    nx.write_graphml(tnet.G_supergraph, "tmp/G_supergraph.txt")
    if exogenous_G != False:
        nx.write_graphml(exogenous_G, "tmp/exogenous_G.txt")
    
    f = open("tmp/g.txt","w")
    f.write( str(tnet.g) )
    f.close()

    f = open("tmp/fcoeffs.txt","w")
    f.write( str(tnet.fcoeffs) )
    f.close()


def solve_social_Juliajson(tnet, exogenous_G=False):
    # Save to json files
    nx2json(tnet.G, "tmp/G.json")
    nx2json(tnet.G_supergraph, "tmp/G_supergraph.json")
    if exogenous_G != False:
        nx2json(exogenous_G, "tmp/exogenous_G.json")

    js = json.dumps({str(k):v for k,v in tnet.g.items()})
    f = open("tmp/g.json", "w")
    f.write(js)
    f.close()

    f = open("tmp/fcoeffs.json", "w")
    f.write(str(tnet.fcoeffs))
    f.close()


'''
import cvxpy as cp
def solve_social_NLP(tnet, exogenous_G=False):

    # Build variables
    xc = {}
    for i,j in tnet.G_supergraph.edges():
        xc[(i,j)] = cp.Variable(name='xc('+str(i) + ','+ str(j) + ")")

    # objective
    if exogenous_G != False:
        obj = 0
        for i,j in tnet.G_supergraph.edges():
            for n in range(len(tnet.fcoeffs)):
                obj += tnet.G_supergraph[i][j]['t_0']*xc[(i,j)]*tnet.fcoeffs[n]*cp.power(xc[(i,j)]+exogenous_G[i][j], n)
    else:
        obj = 0
        for i,j in tnet.G_supergraph.edges():
            for n in range(len(tnet.fcoeffs)):
                obj += tnet.G_supergraph[i][j]['t_0']*xc[(i,j)]*tnet.fcoeffs[n]*cp.power(xc[(i,j)], n)

    cp.Minimize(obj)
    # constraints


    cp.Problem(cp.Minimize(obj)).solve(verbose=True)

    print(xc.values())

'''