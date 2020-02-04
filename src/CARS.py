
from gurobipy import *
import pwlf as pw
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

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

    #plot the results
    #plt.figure()
    #plt.plot(x, y, '-')
    #plt.plot(xHat, yHat, '-')
    #plt.show()
    return  my_pw


def solve_CARS(tnet, exogenous_G, fcoeffs):

    fun = get_approx_fun(fcoeffs)
    beta = get_beta(fun)
    theta = get_theta(fun)

    m = Model('QP')
    m.setParam('OutputFlag',0)


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
            obj += tnet.G_supergraph[i][j]['t_0'] * m.getVarByName('x^'+str(w)+'_'+str(i)+'_'+str(j))
    for i,j in tnet.G.edges():
        obj +=  (tnet.G[i][j]['t_0'] * beta[0]/tnet.G[i][j]['capacity']) * m.getVarByName('e^1_'+str(i)+'_'+str(j)) \
            * (m.getVarByName('e^1_'+str(i)+'_'+str(j)) + theta[0] - exogenous_G[i][j]['flow'] +  m.getVarByName('e^2_'+str(i)+'_'+str(j))) \
            + (tnet.G[i][j]['t_0'] * beta[1]/tnet.G[i][j]['capacity']) * m.getVarByName('e^2_'+str(i)+'_'+str(j)) \
            * (m.getVarByName('e^1_'+str(i)+'_'+str(j)) + theta[1] - exogenous_G[i][j]['flow'])

    m.update()

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


    for j in tnet.G.nodes():
        m.addConstr(quicksum(m.getVarByName('x^R'+str(i)+'_'+str(j)) + quicksum(m.getVarByName('x^' + str(w) + '_' + str(i) + '_' + str(j)) for w in tnet.g.keys() ) for i,l in tnet.G_supergraph.in_edges(nbunch=j)) \
                    == quicksum(m.getVarByName('x^R'+str(j)+'_'+str(k)) + quicksum(m.getVarByName('x^' + str(w) + '_' + str(j) + '_' + str(k)) for w in tnet.g.keys() ) for l,k in tnet.G_supergraph.out_edges(nbunch=j)) )
    m.update()


    for i,j in tnet.G.edges():
        m.addConstr(m.getVarByName('e^1_'+str(i)+'_'+str(j))  >= quicksum(m.getVarByName('x^' + str(w) + '_' + str(i) + '_' + str(j)) for w in tnet.g.keys()) +  m.getVarByName('x^R'+str(i)+'_'+str(j)) + exogenous_G[i][j]['flow'] - theta[0] - m.getVarByName('e^2_'+str(i)+'_'+str(j)))
        m.addConstr(m.getVarByName('e^2_'+str(i)+'_'+str(j)) >= quicksum(m.getVarByName('x^' + str(w) + '_' + str(i) + '_' + str(j)) for w in tnet.g.keys()) +  m.getVarByName('x^R'+str(i)+'_'+str(j)) + exogenous_G[i][j]['flow']- theta[1])

    m.setObjective(obj, GRB.MINIMIZE)
    m.update()
    #m.write("out.lp")
    print(m)
    m.optimize()

    print(obj.getValue())

    for v in m.getVars():
        if v.X>0:
            print('%s = %g' % (v.varName, v.X))


def solve_CARS2(tnet, exogenous_G, fcoeffs, xa=0.01):

    fun = get_approx_fun(fcoeffs, xa=xa, nlines=2)
    beta = get_beta(fun)
    theta = get_theta(fun)
    
    # Set obj func parameters
    Vt = 100
    Vd = 0
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

    # Set exogenous flow
    if exogenous_G == 0:
        exogenous_G = tnet.G.copy()
        for i,j in tnet.G.edges():
            exogenous_G[i][j]['flow'] = 0

    m = Model('QP')
    m.setParam('OutputFlag',0)

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
            * (m.getVarByName('e^1_'+str(i)+'_'+str(j)) + theta[0] - exogenous_G[i][j]['flow']) \
            + (Vd * tnet.G[i][j]['length'] + Ve * tnet.G[i][j]['e']) * (sum(m.getVarByName('x^'+str(w)+'_'+str(i)+'_'+str(j)) for w in tnet.g.keys()) + m.getVarByName('x^R'+str(i)+'_'+str(j)))
    m.update()

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


    for j in tnet.G.nodes():
        m.addConstr(quicksum(m.getVarByName('x^R'+str(i)+'_'+str(j)) + quicksum(m.getVarByName('x^' + str(w) + '_' + str(i) + '_' + str(j)) for w in tnet.g.keys() ) for i,l in tnet.G_supergraph.in_edges(nbunch=j)) \
                    == quicksum(m.getVarByName('x^R'+str(j)+'_'+str(k)) + quicksum(m.getVarByName('x^' + str(w) + '_' + str(j) + '_' + str(k)) for w in tnet.g.keys() ) for l,k in tnet.G_supergraph.out_edges(nbunch=j)) )
    m.update()


    for i,j in tnet.G.edges():
        m.addConstr(m.getVarByName('e^1_'+str(i)+'_'+str(j))  >= quicksum(m.getVarByName('x^' + str(w) + '_' + str(i) + '_' + str(j)) for w in tnet.g.keys()) +  m.getVarByName('x^R'+str(i)+'_'+str(j)) + exogenous_G[i][j]['flow'] - theta[0])

    m.setObjective(obj, GRB.MINIMIZE)
    m.update()
    #m.write("out.lp")
    m.optimize()

    # Summariazing results
    for i,j in tnet.G_supergraph.edges():
        if isinstance(i,str) or isinstance(j,str):
            tnet.G_supergraph[i][j]['flow'] = sum(m.getVarByName('x^' + str(w) + '_' + str(i) + '_' + str(j)).X for w in tnet.g.keys())
            tnet.G_supergraph[i][j]['flowNoRebalancing'] = sum(m.getVarByName('x^' + str(w) + '_' + str(i) + '_' + str(j)).X for w in tnet.g.keys())
        else:
            tnet.G_supergraph[i][j]['flow'] = m.getVarByName('x^R' + str(i) + '_' + str(j)).X + sum(m.getVarByName('x^' + str(w) + '_' + str(i) + '_' + str(j)).X for w in tnet.g.keys())
            tnet.G_supergraph[i][j]['flowNoRebalancing'] = sum(m.getVarByName('x^' + str(w) + '_' + str(i) + '_' + str(j)).X for w in tnet.g.keys())
    tnet.cars_obj = obj.getValue()
    return tnet




def solve_CARS2_noRebalancing(tnet, exogenous_G, fcoeffs, xa=1.2):

    fun = get_approx_fun(fcoeffs, xa=xa, nlines=2)
    beta = get_beta(fun)
    theta = get_theta(fun)

    if exogenous_G == 0:
        exogenous_G = tnet.G.copy()
        for i,j in tnet.G.edges():
            exogenous_G[i][j]['flow'] = 0

    m = Model('QP')
    m.setParam('OutputFlag',0)

    # Define variables
    var_xw = [m.addVar(lb=0, name='x^'+str(w)+'_'+str(i)+'_'+str(j)) for i,j in tnet.G_supergraph.edges() for w in tnet.g.keys()]
    var_e1 = [m.addVar(lb=0, name='e^1_'+str(i)+'_'+str(j)) for i,j in tnet.G.edges()]

    m.update()

    # Set objective
    obj = 0
    for i,j in tnet.G_supergraph.edges():
        for w in tnet.g.keys():
            obj += tnet.G_supergraph[i][j]['t_0'] * m.getVarByName('x^'+str(w)+'_'+str(i)+'_'+str(j))
    for i,j in tnet.G.edges():
        obj +=  (tnet.G[i][j]['t_0'] * beta[0]/tnet.G[i][j]['capacity']) * m.getVarByName('e^1_'+str(i)+'_'+str(j)) \
            * (m.getVarByName('e^1_'+str(i)+'_'+str(j)) + theta[0] - exogenous_G[i][j]['flow'])
    m.update()

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

    for i,j in tnet.G.edges():
        m.addConstr(m.getVarByName('e^1_'+str(i)+'_'+str(j))  >= quicksum(m.getVarByName('x^' + str(w) + '_' + str(i) + '_' + str(j)) for w in tnet.g.keys()) + exogenous_G[i][j]['flow'] - theta[0])

    m.setObjective(obj, GRB.MINIMIZE)
    m.update()
    m.write("out.lp")
    m.optimize()

    # Summariazing results
    for i,j in tnet.G_supergraph.edges():
        if isinstance(i,str) or isinstance(j,str):
            tnet.G_supergraph[i][j]['flow'] = sum(m.getVarByName('x^' + str(w) + '_' + str(i) + '_' + str(j)).X for w in tnet.g.keys())
            tnet.G_supergraph[i][j]['flowNoRebalancing'] = tnet.G_supergraph[i][j]['flow']
        else:
            tnet.G_supergraph[i][j]['flow'] = sum(m.getVarByName('x^' + str(w) + '_' + str(i) + '_' + str(j)).X for w in tnet.g.keys())
            tnet.G_supergraph[i][j]['flowNoRebalancing'] = tnet.G_supergraph[i][j]['flow']
    return tnet


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


def travel_time(tnet, i, j):
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
        [tnet.fcoeffs[n] * (tnet.G_supergraph[i][j]['flow'] / tnet.G_supergraph[i][j]['capacity']) ** n for n in range(len(tnet.fcoeffs))])


def get_totalTravelTime(tnet):
    """
    evalute the travel time function on the SuperGraph level

    Parameters
    ----------

    tnet: transportation network object

    Returns
    -------
    float

    """
    return sum([tnet.G_supergraph[i][j]['flow'] * tnet.G_supergraph[i][j]['t_0'] * travel_time(tnet, i, j) for i, j in tnet.G_supergraph.edges()])


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

def get_totalTravelTime_without_Rebalancing(tnet, G_exogenous):
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


def add_G_flows_no_rebalancing(array):
    # TODO: add description

    G = array[0].copy()
    for tn in array[1:]:
        for i, j in G.edges():
            G[i][j]['flow'] += tn[i][j]['flowNoRebalancing']
    return G