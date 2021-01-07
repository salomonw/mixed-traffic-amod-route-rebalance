
from gurobipy import *
#import pwlf as pw
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from networkx.readwrite import json_graph
import json
from src.utils import *
import src.pwapprox as pw


def get_theta(fun):
    return fun.fit_breaks[1:-1]

def get_beta(fun):
    return fun.slopes[1:]


def eval_travel_time(x, fcoeffs):
    return sum([fcoeffs[i]*x**i for i in range(len(fcoeffs))])


def eval_pw(a, b, theta, x):
    theta2  = theta.copy()
    theta2.append(1000)
    for i in range(len(theta2)-1):
        if (theta2[i]<=x) and (theta2[i+1]>x):
            y =  b[i]+a[i]*x
    return  y



def get_approx_fun(fcoeffs, range_=[0,2], nlines=3, theta=False, plot=False):
    # Generate data
    x = [i  for i in list(np.linspace(range_[0], range_[1], 100))]
    y = [eval_travel_time(i, fcoeffs) for i in x]
    if theta==False:
        pws = pw.pwapprox(x, y, k=nlines)
        pws.fit_convex_boyd(N=30, L=30)
        rms = min(pws.rms_vec)
        i = pws.rms_vec.index(rms)
        a = pws.a_list[i]
        b = pws.b_list[i]
        theta = pws.thetas[i]
        theta.insert(0,0)
        theta.append(range_[1])
    else:
        pws = pw.pwapprox(x, y, k=nlines)
        pws.fit_convex_with_theta(theta)
        theta = theta
        a = pws.a
        rms = 0

    if plot == True:
        fig, ax = plt.subplots(2)
        ax[0].plot(x, y , label = 'Original', color='k')
        ypws = [eval_pw(a,b, theta[0:-1], i) for i in x]
        ax[0].plot(x, ypws, label='pwlinear', color='red')
        for th in theta:
            ax[0].axvline(x=th, linestyle=':')
        plt.grid()
        plt.xlabel('x')
        plt.ylabel('t(x)')
        plt.legend()
        plt.tight_layout()
        pws.plot_rms(ax=ax[1])
        plt.show()
    return  theta, a, rms


@timeit
def add_demand_cnstr(m, tnet, x, bush=False):
    # Set Constraints
    if bush==False:
        for j in tnet.G_supergraph.nodes():
            for w, d in tnet.g.items():
                if j == w[0]:
                    m.addConstr(quicksum(m.getVarByName('x^'+str(w)+'_'+str(i)+'_'+str(j)) for i,l in tnet.G_supergraph.in_edges(nbunch=j)) + d == quicksum(m.getVarByName('x^'+str(w)+'_'+str(j)+'_'+str(k)) for l,k in tnet.G_supergraph.out_edges(nbunch=j)))
                elif j == w[1]:
                    m.addConstr(quicksum(m.getVarByName('x^' + str(w) + '_' + str(i) + '_' + str(j)) for i,l in tnet.G_supergraph.in_edges(nbunch=j)) == quicksum(m.getVarByName('x^' + str(w) + '_' + str(j) + '_' + str(k)) for l,k in tnet.G_supergraph.out_edges(nbunch=j)) + d)
                else:
                    m.addConstr(quicksum(m.getVarByName('x^' + str(w) + '_' + str(i) + '_' + str(j)) for i,l in tnet.G_supergraph.in_edges(nbunch=j)) == quicksum(m.getVarByName('x^' + str(w) + '_' + str(j) + '_' + str(k)) for l,k in tnet.G_supergraph.out_edges(nbunch=j)))
    else:
        p = {j:0  for j in tnet.G_supergraph.nodes()}
        for j in tnet.O:
            p[j] = sum([tnet.g[(s,t)] for s,t in tnet.g.keys() if t==j]) - sum([tnet.g[(s,t)] for s,t in tnet.g.keys() if s==j])

        # Global
        #[m.addConstr(quicksum([x[(i,j)] for i,l in tnet.G_supergraph.in_edges(nbunch=j)]) - quicksum([x[(j,k)] for l,k in tnet.G_supergraph.out_edges(nbunch=j)]) == p[j] ) for j in tnet.G_supergraph.nodes()]
        # Local
        #'''
        dsum = {s:sum([v for k,v in tnet.g.items() if k[0]==s]) for s in tnet.O}
        D = {s:list([d for o, d in tnet.g.keys() if o == s]) for s in tnet.O}
        #l = {s: [j for j in tnet.G_supergraph.nodes() if j not in set(D[s]) if j != s] for s in tnet.O}
        [m.addConstr(quicksum(m.getVarByName('x^' + str(s) + '_' + str(i) + '_' + str(j)) for i, l in tnet.G_supergraph.in_edges(nbunch=j)) - tnet.g[(s,j)] \
                        == quicksum(m.getVarByName('x^' + str(s) + '_' + str(j) + '_' + str(k)) for l, k in tnet.G_supergraph.out_edges(nbunch=j))) for s in tnet.O for j in D[s]]

        [m.addConstr(quicksum(m.getVarByName('x^' + str(s) + '_' + str(i) + '_' + str(s)) for i, l in tnet.G_supergraph.in_edges(nbunch=s)) \
                        == quicksum(m.getVarByName('x^' + str(s) + '_' + str(s) + '_' + str(k)) for l, k in tnet.G_supergraph.out_edges(nbunch=s)) - dsum[s]) for s in tnet.O]

        [m.addConstr(quicksum(m.getVarByName('x^' + str(s) + '_' + str(i) + '_' + str(j)) for i, l in tnet.G_supergraph.in_edges(nbunch=j)) \
                        == quicksum(m.getVarByName('x^' + str(s) + '_' + str(j) + '_' + str(k)) for l, k in tnet.G_supergraph.out_edges(nbunch=j))) \
                        for s in tnet.O for j in [j for j in tnet.G_supergraph.nodes() if j not in set(D[s]) if j != s]]
        '''
        for s in tnet.O:
            dsum = sum([v for k,v in tnet.g.items() if k[0]==s])
            D = [d for o, d in tnet.g.keys() if o == s]

            for j in tnet.G_supergraph.nodes():    
                if j in D:
                    d1 = tnet.g[(s,j)]
                    m.addConstr(quicksum(m.getVarByName('x^' + str(s) + '_' + str(i) + '_' + str(j)) for i, l in tnet.G_supergraph.in_edges(nbunch=j)) - d1 \
                        == quicksum(m.getVarByName('x^' + str(s) + '_' + str(j) + '_' + str(k)) for l, k in tnet.G_supergraph.out_edges(nbunch=j)))
                elif j == s:
                    m.addConstr(quicksum(m.getVarByName('x^' + str(s) + '_' + str(i) + '_' + str(j)) for i, l in tnet.G_supergraph.in_edges(nbunch=j)) \
                        == quicksum(m.getVarByName('x^' + str(s) + '_' + str(j) + '_' + str(k)) for l, k in tnet.G_supergraph.out_edges(nbunch=j)) - dsum)
                else:
                    m.addConstr(quicksum(m.getVarByName('x^' + str(s) + '_' + str(i) + '_' + str(j)) for i, l in tnet.G_supergraph.in_edges(nbunch=j)) \
                                == quicksum(m.getVarByName('x^' + str(s) + '_' + str(j) + '_' + str(k)) for l, k in tnet.G_supergraph.out_edges(nbunch=j)) )
        '''
    m.update()

@timeit
def add_rebalancing_cnstr(m, tnet, xu):
    [m.addConstr(quicksum(m.getVarByName('x^R' + str(i) + '_' + str(j)) + xu[(i, j)] for i, l in tnet.G.in_edges(nbunch=j)) \
        == quicksum(m.getVarByName('x^R' + str(j) + '_' + str(k)) + xu[j, k] for l, k in tnet.G.out_edges(nbunch=j))) for j in tnet.G.nodes()]

    #[m.addConstr(m.getVarByName('x^R'+str(i)+'_'+str(j))==0) for i,j in tnet.G_supergraph.edges() if (type(i)!=int) or (type(j)!=int)]

    m.update()

@timeit
def set_optimal_flows(m , tnet, G_exogenous=False, bush=False):
    if bush:
        for i,j in tnet.G_supergraph.edges():
            tnet.G_supergraph[i][j]['flowNoRebalancing'] = sum(m.getVarByName('x^' + str(s) + '_' + str(i) + '_' + str(j)).X for s in tnet.O)
            tnet.G_supergraph[i][j]['flow'] = tnet.G_supergraph[i][j]['flowNoRebalancing']
            if isinstance(i, int) and isinstance(j, int):
                tnet.G_supergraph[i][j]['flowRebalancing'] = m.getVarByName('x^R' + str(i) + '_' + str(j)).X
                tnet.G_supergraph[i][j]['flow'] += tnet.G_supergraph[i][j]['flowRebalancing']
            #else:
            #tnet.G_supergraph[i][j]['flow'] = tnet.G_supergraph[i][j]['flowRebalancing'] + tnet.G_supergraph[i][j]['flowNoRebalancing']
            tnet.G_supergraph[i][j]['t_k'] = travel_time(tnet, i, j, G_exo=G_exogenous)
    else:
        for i,j in tnet.G_supergraph.edges():
            tnet.G_supergraph[i][j]['flowNoRebalancing'] = sum(m.getVarByName('x^' + str(w) + '_' + str(i) + '_' + str(j)).X for w in tnet.g.keys())
            tnet.G_supergraph[i][j]['flow'] = tnet.G_supergraph[i][j]['flowNoRebalancing']
            if isinstance(i, int) and isinstance(j, int):
                tnet.G_supergraph[i][j]['flowRebalancing'] = m.getVarByName('x^R' + str(i) + '_' + str(j)).X
                tnet.G_supergraph[i][j]['flow'] += tnet.G_supergraph[i][j]['flowRebalancing']

            #else:
            #    tnet.G_supergraph[i][j]['flow'] = m.getVarByName('x^R' + str(i) + '_' + str(j)).X + sum(m.getVarByName('x^' + str(w) + '_' + str(i) + '_' + str(j)).X for w in tnet.g.keys())
            #tnet.G_supergraph[i][j]['flowRebalancing'] = m.getVarByName('x^R' + str(i) + '_' + str(j)).X
            #tnet.G_supergraph[i][j]['flowNoRebalancing'] = sum(m.getVarByName('x^' + str(w) + '_' + str(i) + '_' + str(j)).X for w in tnet.g.keys())
            tnet.G_supergraph[i][j]['t_k'] = travel_time(tnet, i, j, G_exo=G_exogenous)



@timeit
def set_optimal_rebalancing_flows(m,tnet):
    for i,j in tnet.G.edges():
        tnet.G_supergraph[i][j]['flow'] += m.getVarByName('x^R' + str(i) + '_' + str(j)).X
        tnet.G_supergraph[i][j]['flowRebalancing'] = m.getVarByName('x^R' + str(i) + '_' + str(j)).X
        tnet.G_supergraph[i][j]['flowNoRebalancing'] = tnet.G[i][j]['flow'] - tnet.G[i][j]['flowRebalancing']


def eval_obj_funct(tnet, G_exogenous):
    Vt, Vd, Ve = set_CARS_par(tnet)
    obj = Vt * get_totalTravelTime_without_Rebalancing(tnet, G_exogenous=G_exogenous)
    obj = obj + sum([(Vd*tnet.G_supergraph[i][j]['t_0'] + Ve *tnet.G_supergraph[i][j]['e']) * (tnet.G_supergraph[i][j]['flow']-tnet.G_supergraph[i][j]['flowNoRebalancing']) \
               for i,j in tnet.G.edges()])
    return obj/tnet.totalDemand

def set_CARS_par(tnet):
    # Set obj func parameters
    Vt = 24.4
    Vd = 0.286
    Ve = 0.247
    # Set the electricity constant
    ro = 1.25
    Af = 0.4
    cd = 1
    cr = 0.008
    mv = 750
    g = 9.81
    nu = 0.72

    for i,j in tnet.G_supergraph.edges():
        tnet.G_supergraph[i][j]['e'] =  (ro/2 *Af*cd * (tnet.G_supergraph[i][j]['t_0']/tnet.G_supergraph[i][j]['length'])**2 *cr * mv * g)* tnet.G_supergraph[i][j]['length']/nu
    return Vt, Vd, Ve

@timeit
def set_exogenous_flow(tnet, exogenous_G):
    # Set exogenous flow
    exo_G = tnet.G_supergraph.copy()
    for i, j in tnet.G_supergraph.edges():
        exo_G[i][j]['flow'] = 0
    if exogenous_G != False:
        for i,j in exogenous_G.edges():
            exo_G[i][j]['flow'] = exogenous_G[i][j]['flow']
    return  exo_G



'''
    for i,j in tnet.G_supergraph.edges():
        obj += Vt *tnet.G_supergraph[i][j]['t_0'] * quicksum(m.getVarByName('x^'+str(w)+'_'+str(i)+'_'+str(j)) for w in tnet.g.keys()) \
            + Vt * (tnet.G_supergraph[i][j]['t_0'] * beta[0]/tnet.G_supergraph[i][j]['capacity']) * m.getVarByName('e^1_'+str(i)+'_'+str(j)) \
                * (m.getVarByName('e^1_'+str(i)+'_'+str(j)) + theta[0]*tnet.G_supergraph[i][j]['capacity'] - exogenous_G[i][j]['flow']) \
            + Vt * (tnet.G_supergraph[i][j]['t_0'] * beta[1]/tnet.G_supergraph[i][j]['capacity']) * m.getVarByName('e^2_'+str(i)+'_'+str(j)) \
                * (m.getVarByName('e^2_'+str(i)+'_'+str(j)) + theta[1]*tnet.G_supergraph[i][j]['capacity'] - exogenous_G[i][j]['flow']) \
            + Vt * (tnet.G_supergraph[i][j]['t_0'] * beta[0]/tnet.G_supergraph[i][j]['capacity'] * m.getVarByName('e^2_'+str(i)+'_'+str(j)) \
                * (theta[1]*tnet.G_supergraph[i][j]['capacity'] - theta[0]*tnet.G_supergraph[i][j]['capacity']) )

    for i,j in tnet.G.edges():
        obj +=  (Vd * tnet.G_supergraph[i][j]['t_0'] + Ve * tnet.G_supergraph[i][j]['e']) * ( \
                #sum(m.getVarByName('x^' + str(w) + '_' + str(i) + '_' + str(j)) for w in tnet.g.keys()) +\
                m.getVarByName('x^R' + str(i) + '_' + str(j)))
'''
@timeit
def get_obj_CARSn(m, tnet, xu,  theta, a, exogenous_G, linear=False):#, userCentric=False):
    #TODO: this could be written more efficiently, include user-centric approach
    #if linear:
    #userCentric = False
    #if userCentric != True:
    Vt, Vd, Ve = set_CARS_par(tnet)
    if linear == True:
        obj = quicksum(Vt * tnet.G_supergraph[i][j]['t_0'] * xu[(i, j)] for i,j in tnet.G_supergraph.edges())
        obj += quicksum(\
                quicksum( Vt * tnet.G_supergraph[i][j]['t_0'] * a[l]/tnet.G_supergraph[i][j]['capacity'] *( \
                m.getVarByName('e^'+str(l)+'_'+str(i)+'_'+str(j)) * (0+quicksum(((theta[k + 1] - theta[k])*tnet.G_supergraph[i][j]['capacity']) for k in range(0,l))) \
                + m.getVarByName('e^'+str(l)+'_'+str(i)+'_'+str(j)) * ( (theta[l + 1] - theta[l])*tnet.G_supergraph[i][j]['capacity'] ) \
                + (theta[l+1] - theta[l])*tnet.G_supergraph[i][j]['capacity']*(0+quicksum(m.getVarByName('e^'+str(k)+'_'+str(i)+'_'+str(j)) for k in range(l+1, len(theta)-1))) \
                - m.getVarByName('e^'+str(l)+'_'+str(i)+'_'+str(j)) * exogenous_G[i][j]['flow'] \
                ) for l in range(len(theta)-1))  \
                + (Vd * tnet.G_supergraph[i][j]['t_0'] + Ve * tnet.G_supergraph[i][j]['e']) * m.getVarByName('x^R' + str(i) + '_' + str(j))\
                for i,j in tnet.G.edges())
        '''
        obj = quicksum( Vt * tnet.G_supergraph[i][j]['t_0'] * xu[(i, j)] \
                + quicksum( Vt * tnet.G_supergraph[i][j]['t_0'] * a[l]/tnet.G_supergraph[i][j]['capacity'] *  m.getVarByName('e^'+str(l)+'_'+str(i)+'_'+str(j)) * (quicksum(((theta[k + 1] - theta[k]) * tnet.G_supergraph[i][j]['capacity']) for k in range(l))) \
                + Vt * tnet.G_supergraph[i][j]['t_0'] * a[l]/tnet.G_supergraph[i][j]['capacity'] * m.getVarByName('e^'+str(l)+'_'+str(i)+'_'+str(j)) * ((theta[l+1]-theta[l])*tnet.G_supergraph[i][j]['capacity']) \
                + Vt * tnet.G_supergraph[i][j]['t_0'] * a[l]/tnet.G_supergraph[i][j]['capacity'] *(theta[l+1] - theta[l])*tnet.G_supergraph[i][j]['capacity']*(quicksum(m.getVarByName('e^'+str(k)+'_'+str(i)+'_'+str(j)) for k in range(l+1, len(theta)-1))) \
                - Vt * tnet.G_supergraph[i][j]['t_0'] * a[l]/tnet.G_supergraph[i][j]['capacity'] * m.getVarByName('e^'+str(l)+'_'+str(i)+'_'+str(j)) * exogenous_G[i][j]['flow'] \
                for l in range(len(theta)-1))  \
                + (Vd * tnet.G_supergraph[i][j]['t_0'] + Ve * tnet.G_supergraph[i][j]['e']) * m.getVarByName('x^R' + str(i) + '_' + str(j)) \
                for i,j in tnet.G_supergraph.edges())
        '''
    else:
        obj = quicksum(Vt * tnet.G_supergraph[i][j]['t_0'] * xu[(i, j)] for i,j in tnet.G_supergraph.edges())
        obj = obj+ quicksum(\
                quicksum( Vt * tnet.G_supergraph[i][j]['t_0'] * a[l]/tnet.G_supergraph[i][j]['capacity'] *  m.getVarByName('e^'+str(l)+'_'+str(i)+'_'+str(j)) * (quicksum(((theta[k + 1] - theta[k])*tnet.G_supergraph[i][j]['capacity']) for k in range(0,l))) \
                + Vt * tnet.G_supergraph[i][j]['t_0'] * a[l]/tnet.G_supergraph[i][j]['capacity'] * m.getVarByName('e^'+str(l)+'_'+str(i)+'_'+str(j)) * ( m.getVarByName('e^'+str(l)+'_'+str(i)+'_'+str(j))) \
                + Vt * tnet.G_supergraph[i][j]['t_0'] * a[l]/tnet.G_supergraph[i][j]['capacity'] *(theta[l+1] - theta[l])*tnet.G_supergraph[i][j]['capacity']*(quicksum(m.getVarByName('e^'+str(k)+'_'+str(i)+'_'+str(j)) for k in range(l+1, len(theta)-1))) \
                - Vt * tnet.G_supergraph[i][j]['t_0'] * a[l]/tnet.G_supergraph[i][j]['capacity'] * m.getVarByName('e^'+str(l)+'_'+str(i)+'_'+str(j)) * exogenous_G[i][j]['flow'] \
                for l in range(len(theta)-1)) \
                + (Vd * tnet.G_supergraph[i][j]['t_0'] + Ve * tnet.G_supergraph[i][j]['e']) * m.getVarByName('x^R' + str(i) + '_' + str(j))\
                for i,j in tnet.G.edges())


        '''
        if linear == True:
            Vt, Vd, Ve = set_CARS_par(tnet)
            obj = quicksum( Vt * tnet.G_supergraph[i][j]['t_0'] * xu[(i, j)] \
                    + quicksum( Vt * tnet.G_supergraph[i][j]['t_0'] * a[l]/tnet.G_supergraph[i][j]['capacity'] *  m.getVarByName('e^'+str(l)+'_'+str(i)+'_'+str(j)) * (quicksum(((theta[k + 1] - theta[k]) * tnet.G_supergraph[i][j]['capacity']) for k in range(l))) \
                    + Vt * tnet.G_supergraph[i][j]['t_0'] * a[l]/tnet.G_supergraph[i][j]['capacity'] * m.getVarByName('e^'+str(l)+'_'+str(i)+'_'+str(j)) * ( m.getVarByName('e^'+str(l)+'_'+str(i)+'_'+str(j))) \
                    + Vt * tnet.G_supergraph[i][j]['t_0'] * a[l]/tnet.G_supergraph[i][j]['capacity'] *(theta[l+1] - theta[l])*tnet.G_supergraph[i][j]['capacity']*(quicksum(m.getVarByName('e^'+str(k)+'_'+str(i)+'_'+str(j)) for k in range(l+1, len(theta)-1))) \
                    - Vt * tnet.G_supergraph[i][j]['t_0'] * a[l]/tnet.G_supergraph[i][j]['capacity'] * m.getVarByName('e^'+str(l)+'_'+str(i)+'_'+str(j)) * exogenous_G[i][j]['flow'] \
                    for l in range(len(theta)-1))  \
                    + (Vd * tnet.G_supergraph[i][j]['t_0'] + Ve * tnet.G_supergraph[i][j]['e']) * m.getVarByName('x^R' + str(i) + '_' + str(j)) \
                    for i,j in tnet.G_supergraph.edges())

        else:
            Vt, Vd, Ve = set_CARS_par(tnet)
            obj = quicksum(Vt * tnet.G_supergraph[i][j]['t_0'] * xu[(i,j)] for i,j in tnet.G_supergraph.edges())
            for i,j in tnet.G_supergraph.edges():
                t0 = tnet.G_supergraph[i][j]['t_0']
                mij = tnet.G_supergraph[i][j]['capacity']
                ue = exogenous_G[i][j]['flow']
                for l in range(len(theta)-1):
                    Vtt_0al = Vt * t0 * a[l]/mij
                    e_l = m.getVarByName('e^'+str(l)+'_'+str(i)+'_'+str(j))
                    obj += Vtt_0al *  e_l * (quicksum(((theta[k + 1] - theta[k]) * mij) for k in range(l)))
                    obj += Vtt_0al * e_l * ( e_l)
                    obj += Vtt_0al *(theta[l+1] - theta[l])*mij*(quicksum(m.getVarByName('e^'+str(k)+'_'+str(i)+'_'+str(j)) for k in range(l+1, len(theta)-1)))
                    obj -= Vtt_0al* e_l * ue
        '''


    '''   
    else:
        #if linear == True:
        Vt, Vd, Ve = set_CARS_par(tnet)
        obj = 0
        obj = quicksum(Vt * tnet.G_supergraph[i][j]['t_0'] for i,j in tnet.G_supergraph.edges())
        for i,j in tnet.G_supergraph.edges():
            t0 = tnet.G_supergraph[i][j]['t_0']
            mij = tnet.G_supergraph[i][j]['capacity']
            ue = exogenous_G[i][j]['flow']
            for l in range(len(theta)-1):
                Vtt_0al = Vt * t0 * a[l] / mij
                e_l = m.getVarByName('e^' + str(l) + '_' + str(i) + '_' + str(j))
                obj += quicksum(Vtt_0al * e_l / ((theta[k + 1] - theta[k]) * mij - ue) for k in range(l))
                obj += Vtt_0al * (1 - e_l / ue)                             #Vtt_0al * (e_l * e_l - e_l * ue)
                obj += quicksum(Vtt_0al * ((theta[l + 1] - theta[l]) * mij / m.getVarByName(
                    'e^' + str(k) + '_' + str(i) + '_' + str(j)) - e_l * ue) for k in range(l + 1, len(theta) - 1))
    '''

    #obj += quicksum((Vd * tnet.G_supergraph[i][j]['t_0'] + Ve * tnet.G_supergraph[i][j]['e']) * m.getVarByName('x^R' + str(i) + '_' + str(j)) for i,j in tnet.G_supergraph.edges())
    return obj






@timeit
def add_epsilon_cnstr(m, tnet, xu, n, theta, exogenous_G):
    [m.addConstr(m.getVarByName('e^'+str(l)+'_'+str(i)+'_'+str(j))\
                      >=  xu[(i,j)] \
                      +  m.getVarByName('x^R' + str(i) + '_' + str(j)) \
                      + exogenous_G[i][j]['flow'] \
                      - theta[l]*tnet.G_supergraph[i][j]['capacity'] \
                      - quicksum(m.getVarByName('e^'+str(l+k+1)+'_'+str(i)+'_'+str(j)) for k in range(n-l-1))) for i,j in tnet.G.edges() for l in range(n) ]
                      #- quicksum(m.getVarByName('e^' + str(l + k) + '_' + str(i) + '_' + str(j)) for k in range(n - l))) for i, j in tnet.G_supergraph.edges() for l in range(n)]
    #[m.addConstr(m.getVarByName('e^'+str(l)+'_'+str(i)+'_'+str(j))\
    #                  <= (theta[l+1]-theta[l])*tnet.G_supergraph[i][j]['capacity'] )for i,j in tnet.G_supergraph.edges() for l in range(n-1)]
    # maximum flow
    #[m.addConstr(m.getVarByName('e^' + str(l) + '_' + str(i) + '_' + str(j)) \
    #             <= (theta[l + 1] - theta[l]) * tnet.G_supergraph[i][j]['capacity']) for i, j in
    # tnet.G_supergraph.edges() for l in range(n - 1)]

'''
@timeit
def solve_CARSn(tnet, fcoeffs, n=3, exogenous_G=False, rebalancing=True, linear=False, method=-1, theta=False, a=False):
    if (theta==False) or (a==False):
        systemOptimal = True
        if systemOptimal:
            fcoeffs = SO_fcoeffs(fcoeffs)
            print(fcoeffs)
        theta, a, rms  = get_approx_fun(fcoeffs=fcoeffs, nlines=n, range_=[0,3], plot=False)
        theta.append(3)
    exogenous_G = set_exogenous_flow(tnet, exogenous_G)
    # Start model
    m = Model('QP')
    m.setParam('OutputFlag',0)
    m.setParam('BarHomogeneous', 1)
    m.setParam('Method', method)
    m.update()
    # Define variables
    [m.addVar(lb=0, name='x^'+str(w)+'_'+str(i)+'_'+str(j)) for i,j in tnet.G_supergraph.edges() for w in tnet.g.keys()]
    if rebalancing == True:
        [m.addVar(lb=0, name='x^R'+str(i)+'_'+str(j)) for i,j in tnet.G_supergraph.edges()]
    else:
        [m.addVar(lb=0, ub=0, name='x^R' + str(i) + '_' + str(j)) for i, j in tnet.G_supergraph.edges()]

    [m.addVar(lb=0, name='e^'+str(l)+'_'+str(i)+'_'+str(j)) for i,j in tnet.G_supergraph.edges() for l in range(n)]
    m.update()
    # Set Obj
    obj = get_obj_CARSn(m, tnet, theta, a, exogenous_G, linear=linear)
    # Set Constraints
    add_demand_cnstr(m, tnet)
    if rebalancing==True:
        add_rebalancing_cnstr(m,tnet)
    add_epsilon_cnstr(m, tnet, n, theta, exogenous_G)
    m.update()
    # Solve problem
    m.setObjective(obj, GRB.MINIMIZE)
    m.update()
    m.optimize()
    # saving  results
    set_optimal_flows(m, tnet)
    tnet.cars_obj = obj.getValue()
    od_flows = get_OD_result_flows(m, tnet)
    return tnet, m.Runtime, od_flows
'''

@timeit
def solve_bush_CARSn(tnet, fcoeffs, n=3, exogenous_G=False, rebalancing=True, linear=False, LP_method=-1, QP_method=-1, theta=False, a=False, bush=False, theta_n=3, userCentric=False, od_flows_flag=True):
    #TODO: implement option to select between origin or destination
    fc = fcoeffs.copy()
    if (theta==False) or (a==False):
        if userCentric:
            #fc=fc
            #fc.insert(0,0)
            fc = UC_fcoeffs(fc)
            #print(fc)
        #else:
            #fc.insert(0, 0)
        theta, a, rms  = get_approx_fun(fcoeffs=fc, nlines=n, range_=[0,theta_n], plot=False)
        #a.append(a[-1])
    exogenous_G = set_exogenous_flow(tnet, exogenous_G)
    # Start model
    m = Model('CARS')
    m.setParam('OutputFlag',0 )
    m.setParam('BarHomogeneous', 1)
    #m.setParam("LogToConsole", 0)
    #m.setParam("CSClientLog", 0)
    if linear:
        m.setParam('Method', LP_method)
    else:
        m.setParam('Method', QP_method)
    m.update()

    # Find origins
    tnet.O = list(set([w[0] for w, d in tnet.g.items() if d > 0]))

    # Define variables
    if bush == True:
        [m.addVar(lb=0, name='x^'+str(s)+'_'+str(i)+'_'+str(j)) for i,j in tnet.G_supergraph.edges() for s in tnet.O]
    else:
        [m.addVar(lb=0, name='x^' + str(w) + '_' + str(i) + '_' + str(j)) for i, j in tnet.G_supergraph.edges() for w, d in tnet.g.items()]

    m.update()

    if userCentric==True:
        for i, j in tnet.G_supergraph.edges():
            if isinstance(i, int) and isinstance(j, int):
                continue
            else:
                for s in tnet.O:
                    m.addConstr(m.getVarByName('x^'+str(s)+'_'+str(i)+'_'+str(j)) == 0)

    if rebalancing == True:
        [m.addVar(lb=0, name='x^R'+str(i)+'_'+str(j)) for i,j in tnet.G.edges()]
    else:
        [m.addVar(lb=0, ub=0, name='x^R' + str(i) + '_' + str(j)) for i, j in tnet.G.edges()]

    [m.addVar(name='e^'+str(l)+'_'+str(i)+'_'+str(j), \
              lb=0)# ub=theta[l+1]-theta[l]) \
               for i,j in tnet.G.edges() for l in range(n)]
    #[m.addVar(name='e^'+str(n+1)+'_'+str(i)+'_'+str(j), lb=0) for i,j in tnet.G.edges()]
    m.update()

    if bush==True:
        xu = {(i, j): quicksum(m.getVarByName('x^'+str(s)+'_'+str(i)+'_'+str(j)) for s in tnet.O) for i, j in tnet.G_supergraph.edges()}
    else:
        xu = {(i, j): quicksum(m.getVarByName('x^'+str(w)+'_'+str(i)+'_'+str(j)) for w, d in tnet.g.items()) for i, j in tnet.G_supergraph.edges()}

    # Set Obj
    obj = get_obj_CARSn(m, tnet, xu, theta, a, exogenous_G, linear=linear)

    # Set Constraints
    add_epsilon_cnstr(m, tnet, xu, n, theta, exogenous_G)
    m.update()
    #print(m.display())
    #pause
    add_demand_cnstr(m, tnet, xu,  bush=bush)
    if rebalancing==True:
        add_rebalancing_cnstr(m, tnet, xu)


    # Solve problem
    m.setObjective(obj, GRB.MINIMIZE)
    m.update()
    m.optimize()
    status = {2:'optimal', 3:'infeasible !', 4:'infeasible or unbounded !', 5:'unbounded', 6:'cutoff', 7:'time limit'}
    #print('solver stats: ' + status[GRB.OPTIMAL])

    # saving  results
    set_optimal_flows(m, tnet, G_exogenous=exogenous_G, bush=bush)
    tnet.cars_obj = obj.getValue()
    if od_flows_flag==True:
        od_flows = get_OD_result_flows(m, tnet, bush=bush)
        return tnet, m.Runtime, od_flows
    else:
        return tnet, m.Runtime

'''
def solve_CARS2(tnet, fcoeffs, exogenous_G=False, rebalancing=True):
    #TODO: add description
    fun = get_approx_fun(fcoeffs, range_=[0,2])
    #beta = get_beta(fun)
    #theta = get_theta(fun)
    #theta = [0.9, 1.5] #BPR EMA
    #beta = [1.0, 1.88] #BPR EMA
    #theta = [0.7, 1.2] #BPR NYC
    #beta = [0.5, 1.88] #BPR NYC
    theta = [0.7, 1.2] #estimated NYC
    beta = [0.5, 2.4] #estimated NYC
    #theta = [1.19, 1.2] #CARS
    #beta = [2.44, 2.45] #CARS

    #print(theta)
    #print(beta)

    Vt, Vd, Ve = set_CARS_par(tnet)

    # Set exogenous flow
    if exogenous_G == False:
        exogenous_G = tnet.G_supergraph.copy()
        for i,j in tnet.G_supergraph.edges():
            exogenous_G[i][j]['flow'] = 0

    # Start model
    m = Model('QP')
    m.setParam('OutputFlag',0)
    m.setParam('BarHomogeneous', 1)
    m.setParam('Method', 1)
    m.update()

    # Define variables
    [m.addVar(lb=0, name='x^'+str(w)+'_'+str(i)+'_'+str(j)) for i,j in tnet.G_supergraph.edges() for w in tnet.g.keys()]

    if rebalancing == True:
        [m.addVar(lb=0, name='x^R'+str(i)+'_'+str(j)) for i,j in tnet.G.edges()]
    else:
        [m.addVar(lb=0, ub=0, name='x^R' + str(i) + '_' + str(j)) for i, j in tnet.G.edges()]

    [m.addVar(lb=0, name='e^1_'+str(i)+'_'+str(j)) for i,j in tnet.G_supergraph.edges()]
    [m.addVar(lb=0, name='e^2_'+str(i)+'_'+str(j)) for i,j in tnet.G_supergraph.edges()]

    xu = {(i, j): quicksum(m.getVarByName('x^' + str(w) + '_' + str(i) + '_' + str(j)) for w, d in tnet.g.items()) for i, j in tnet.G_supergraph.edges()}

    m.update()

    # Set objective
    obj = 0

    for i,j in tnet.G_supergraph.edges():
        obj += Vt *tnet.G_supergraph[i][j]['t_0'] * quicksum(m.getVarByName('x^'+str(w)+'_'+str(i)+'_'+str(j)) for w in tnet.g.keys()) \
            + Vt * (tnet.G_supergraph[i][j]['t_0'] * beta[0]/tnet.G_supergraph[i][j]['capacity']) * m.getVarByName('e^1_'+str(i)+'_'+str(j)) \
                * (m.getVarByName('e^1_'+str(i)+'_'+str(j)) + theta[0]*tnet.G_supergraph[i][j]['capacity'] - exogenous_G[i][j]['flow']) \
            + Vt * (tnet.G_supergraph[i][j]['t_0'] * beta[1]/tnet.G_supergraph[i][j]['capacity']) * m.getVarByName('e^2_'+str(i)+'_'+str(j)) \
                * (m.getVarByName('e^2_'+str(i)+'_'+str(j)) + theta[1]*tnet.G_supergraph[i][j]['capacity'] - exogenous_G[i][j]['flow']) \
            + Vt * (tnet.G_supergraph[i][j]['t_0'] * beta[0]/tnet.G_supergraph[i][j]['capacity'] * m.getVarByName('e^2_'+str(i)+'_'+str(j)) \
                * (theta[1]*tnet.G_supergraph[i][j]['capacity'] - theta[0]*tnet.G_supergraph[i][j]['capacity']) )

    for i,j in tnet.G.edges():
        obj +=  (Vd * tnet.G_supergraph[i][j]['t_0'] + Ve * tnet.G_supergraph[i][j]['e']) * ( \
                #sum(m.getVarByName('x^' + str(w) + '_' + str(i) + '_' + str(j)) for w in tnet.g.keys()) +\
                m.getVarByName('x^R' + str(i) + '_' + str(j)))
    m.update()

    # Set Constraints
    add_demand_cnstr(m, tnet, xu)
    if rebalancing==True:
        add_rebalancing_cnstr(m,tnet, xu)
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
    od_flows = get_OD_result_flows(m, tnet)
    return tnet, m.Runtime, od_flows



def solve_CARS(tnet, fcoeffs, exogenous_G=False, rebalancing=True, xa=0.01):
    fun = get_approx_fun(fcoeffs, xa=xa, nlines=2, range_=[0,2])
    beta = get_beta(fun)
    theta = get_theta(fun)

    #print(beta)
    #print(theta)

    Vt, Vd, Ve = set_CARS_par(tnet)
    # Set exogenous flow
    if exogenous_G == False:
        exogenous_G = tnet.G_supergraph.copy()
        for i,j in tnet.G_supergraph.edges():
            exogenous_G[i][j]['flow'] = 0

    m = Model('QP')
    m.setParam('OutputFlag',0)
    m.setParam('BarHomogeneous', 1)
    m.setParam('Method', 1)
    m.update()

    # Define variables
    [m.addVar(lb=0, name='x^'+str(w)+'_'+str(i)+'_'+str(j)) for i,j in tnet.G_supergraph.edges() for w in tnet.g.keys()]

    if rebalancing == True:
        [m.addVar(lb=0, name='x^R'+str(i)+'_'+str(j)) for i,j in tnet.G.edges()]
    else:
        [m.addVar(lb=0, ub=0, name='x^R' + str(i) + '_' + str(j)) for i, j in tnet.G.edges()]

    [m.addVar(lb=0, name='e^1_'+str(i)+'_'+str(j)) for i,j in tnet.G_supergraph.edges()]
    m.update()

    xu = {(i, j): quicksum(m.getVarByName('x^' + str(w) + '_' + str(i) + '_' + str(j)) for w, d in tnet.g.items()) for
          i, j in tnet.G_supergraph.edges()}

    # Set objective
    obj = 0
    for i,j in tnet.G_supergraph.edges():
        for w in tnet.g.keys():
            obj += Vt * tnet.G_supergraph[i][j]['t_0'] * m.getVarByName('x^'+str(w)+'_'+str(i)+'_'+str(j))
        obj += Vt * (tnet.G_supergraph[i][j]['t_0'] * beta[0]/tnet.G_supergraph[i][j]['capacity']) * m.getVarByName('e^1_'+str(i)+'_'+str(j)) \
            * (m.getVarByName('e^1_'+str(i)+'_'+str(j)) + theta[0]*tnet.G_supergraph[i][j]['capacity'] - exogenous_G[i][j]['flow'])
    for i,j in tnet.G.edges():
        obj += (Vd * tnet.G_supergraph[i][j]['t_0'] + Ve * tnet.G_supergraph[i][j]['e']) * ( \
                    #sum(m.getVarByName('x^'+str(w)+'_'+str(i)+'_'+str(j)) for w in tnet.g.keys())  \
                    m.getVarByName('x^R'+str(i)+'_'+str(j)))
    m.update()

    # Set Constraints
    add_demand_cnstr(m, tnet, xu)
    if rebalancing==True:
        add_rebalancing_cnstr(m,tnet, xu)

    for i,j in tnet.G.edges():
        m.addConstr(m.getVarByName('e^1_'+str(i)+'_'+str(j))  >= quicksum(m.getVarByName('x^' + str(w) + '_' + str(i) + '_' + str(j)) for w in tnet.g.keys()) +  m.getVarByName('x^R'+str(i)+'_'+str(j)) + exogenous_G[i][j]['flow'] - theta[0]* tnet.G[i][j]['capacity'])
    m.update()

    m.setObjective(obj, GRB.MINIMIZE)
    m.update()
    m.optimize()
    # saving  results
    set_optimal_flows(m, tnet)
    tnet.cars_obj = obj.getValue()
    return tnet, m.Runtime



def get_obj_CARSn(m, tnet, xu, theta, a, exogenous_G, linear=True, userCentric=False):
    if userCentric != True:
        if linear == True:
            Vt, Vd, Ve = set_CARS_par(tnet)
            obj = quicksum(Vt * tnet.G_supergraph[i][j]['t_0'] * xu[(i,j)] for i,j in tnet.G_supergraph.edges())
            for i,j in tnet.G_supergraph.edges():
                t0 = tnet.G_supergraph[i][j]['t_0']
                for l in range(len(theta)-1):
                    mij = tnet.G_supergraph[i][j]['capacity']
                    ue = exogenous_G[i][j]['flow']
                    Vtt_0al = Vt * t0 * a[l]/mij
                    e_l = m.getVarByName('e^'+str(l)+'_'+str(i)+'_'+str(j))
                    obj += quicksum(Vtt_0al * e_l * ((theta[k + 1] - theta[k]) * mij - ue) for k in range(l))
                    obj += Vtt_0al * (e_l*(theta[l+1]-theta[l])*mij - e_l*ue)
                    obj += quicksum(Vtt_0al * ((theta[l+1] - theta[l])*mij * m.getVarByName('e^'+str(k)+'_'+str(i)+'_'+str(j)) - e_l*ue) for k in range(l+1, len(theta)-1))
        else:
            Vt, Vd, Ve = set_CARS_par(tnet)
            obj = quicksum(Vt * tnet.G_supergraph[i][j]['t_0'] * xu[(i,j)] for i,j in tnet.G_supergraph.edges())
            for i,j in tnet.G_supergraph.edges():
                t0 = tnet.G_supergraph[i][j]['t_0']
                for l in range(len(theta)-1):
                    mij = tnet.G_supergraph[i][j]['capacity']
                    ue = exogenous_G[i][j]['flow']
                    Vtt_0al = Vt * t0 * a[l]/mij
                    e_l = m.getVarByName('e^'+str(l)+'_'+str(i)+'_'+str(j))
                    obj += quicksum(Vtt_0al * e_l * ((theta[k + 1] - theta[k]) * mij - ue) for k in range(l))
                    obj += Vtt_0al * (e_l * e_l - e_l * ue)
                    obj += quicksum(Vtt_0al * ((theta[l+1] - theta[l])*mij * m.getVarByName('e^'+str(k)+'_'+str(i)+'_'+str(j)) - e_l*ue) for k in range(l+1, len(theta)-1))
    else:
        #if linear == True:
        Vt, Vd, Ve = set_CARS_par(tnet)
        obj = 0
        obj = quicksum(Vt * tnet.G_supergraph[i][j]['t_0'] for i,j in tnet.G_supergraph.edges())
        for i,j in tnet.G_supergraph.edges():
            t0 = tnet.G_supergraph[i][j]['t_0']
            mij = tnet.G_supergraph[i][j]['capacity']
            ue = exogenous_G[i][j]['flow']
            for l in range(len(theta)-1):
                Vtt_0al = Vt * t0 * a[l] / mij
                e_l = m.getVarByName('e^' + str(l) + '_' + str(i) + '_' + str(j))
                obj += quicksum(Vtt_0al * e_l / ((theta[k + 1] - theta[k]) * mij - ue) for k in range(l))
                obj += Vtt_0al * (1 - e_l / ue)#Vtt_0al * (e_l * e_l - e_l * ue)
                obj += quicksum(Vtt_0al * ((theta[l + 1] - theta[l]) * mij / m.getVarByName(
                    'e^' + str(k) + '_' + str(i) + '_' + str(j)) - e_l * ue) for k in range(l + 1, len(theta) - 1))


    obj += quicksum((Vd * tnet.G_supergraph[i][j]['t_0'] + Ve * tnet.G_supergraph[i][j]['e']) * m.getVarByName('x^R' + str(i) + '_' + str(j)) for i,j in tnet.G_supergraph.edges())
    return obj
'''


def get_CARS_obj_val(tnet, G_exogenous):
    Vt, Vd, Ve = set_CARS_par(tnet)
    tt = get_totalTravelTime_without_Rebalancing(tnet, G_exogenous=G_exogenous)
    reb = get_rebalancing_total_cost(tnet)
    obj = Vt * tt + reb
    return obj

def get_rebalancing_total_cost(tnet):
    Vt, Vd, Ve = set_CARS_par(tnet)
    reb = get_rebalancing_flow(tnet)
    obj = sum((Vd * tnet.G_supergraph[i][j]['t_0'] + Ve * tnet.G_supergraph[i][j]['e']) * tnet.G_supergraph[i][j]['flowRebalancing']  for i,j in tnet.G_supergraph.edges())
    return obj

@timeit
def get_OD_result_flows(m, tnet, bush=False):
    dic = {}
    if bush==True:
        for s in tnet.O:
            dic[s] = {}
            for i,j in tnet.G_supergraph.edges():
                dic[s][(i,j)] = m.getVarByName('x^' + str(s) + '_' + str(i) + '_' + str(j)).X
    else:
        for w in tnet.g.keys():
            dic[w] = {}
            for i,j in tnet.G_supergraph.edges():
                dic[w][(i,j)] = m.getVarByName('x^' + str(w) + '_' + str(i) + '_' + str(j)).X
    return dic


def solve_rebalancing(tnet, exogenous_G=0):
    Vt, Vd, Ve = set_CARS_par(tnet)
    # Set exogenous flow
    if exogenous_G == 0:
        exogenous_G = tnet.G.copy()
        for i, j in tnet.G.edges():
            exogenous_G[i][j]['flow'] = 0

    m = Model('QP')
    m.setParam('OutputFlag', 1)
    m.setParam('BarHomogeneous', 0)
    m.setParam('Method', 1)
    m.update()

    # Define variables
    [m.addVar(lb=0, name='x^R' + str(i) + '_' + str(j)) for i, j in tnet.G.edges()]
    m.update()

    # Set objective
    #obj = quicksum((Vd * tnet.G[i][j]['t_0'] + Ve * tnet.G[i][j]['e']) * m.getVarByName('x^R' + str(i) + '_' + str(j)) for i,j in tnet.G.edges())
    obj = quicksum((Vt * tnet.G[i][j]['t_k']) * m.getVarByName('x^R' + str(i) + '_' + str(j)) for i, j in
        tnet.G.edges())
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
    return m.Runtime

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
        [tnet.fcoeffs[n] * ((tnet.G_supergraph[i][j]['flowNoRebalancing'] +exo )/ tnet.G_supergraph[i][j]['capacity']) ** n for n in range(len(tnet.fcoeffs))])

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
            't_0'] * travel_time(tnet, i, j) for i, j in
                    tnet.G_supergraph.edges()])
    else:
        return sum([tnet.G_supergraph[i][j]['flowNoRebalancing'] * tnet.G_supergraph[i][j][
            't_0'] * travel_time(tnet, i, j, G_exo=G_exogenous) for i, j in
                    tnet.G_supergraph.edges()])



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
    return sum([tnet.G_supergraph[i][j]['flow']*tnet.G_supergraph[i][j]['length'] for i,j in tnet.G_supergraph.edges() if tnet.G_supergraph[i][j]['type']=='p'])


def get_layer_flow(tnet, symb="'"):
    """
    get  flow in a layer of supergraph

    Parameters
    ----------

    tnet: transportation network object

    Returns
    -------
    float

    """
    return sum([tnet.G_supergraph[i][j]['flow']*tnet.G_supergraph[i][j]['length'] for i,j in tnet.G_supergraph.edges() if tnet.G_supergraph[i][j]['type']==symb])


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
    return sum([tnet.G_supergraph[i][j]['flowNoRebalancing']*tnet.G_supergraph[i][j]['length'] for i,j in tnet.G.edges()])

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
    return sum([(tnet.G_supergraph[i][j]['flow']-tnet.G_supergraph[i][j]['flowNoRebalancing'])*tnet.G_supergraph[i][j]['length'] for i,j in tnet.G.edges()])


def UC_fcoeffs(fcoeffs):
    f = [fcoeffs[i]/(i+1) for i in range(len(fcoeffs))]
    return f

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
    #tnet.G_supergraph = tnet.G
    for i,j in tnet.G_supergraph.edges():
        try:
            tnet.G_supergraph[i][j]['flow'] = tnet.G[i][j]['flow']
        except:
            tnet.G_supergraph[i][j]['flow'] = 0
        tnet.G_supergraph[i][j]['flowNoRebalancing'] = tnet.G_supergraph[i][j]['flow']

def add_G_flows_no_rebalancing(array):
    # TODO: add description

    G = array[0].copy()
    for tn in array[1:]:
        for i, j in G.edges():
            G[i][j]['flow'] += tn[i][j]['flowNoRebalancing']
    return G


def solveMSAsocialCARS(tnet, exogenous_G=False):
    runtime = time.process_time()
    if exogenous_G == False:
        tnet.solveMSAsocial_supergraph()
    else:
        tnet.solveMSAsocial_supergraph(exogenous_G=exogenous_G)
    t = time.process_time() - runtime
    G2supergraph(tnet)
    return t, tnet.TAP.RG


def nx2json(G, fname, exo=False):
    if exo==False:
        D = G.copy()
        for i,j in D.edges():
            D[i][j]['flow'] = 0
        with open(fname, 'w') as outfile1:
            outfile1.write(json.dumps(json_graph.node_link_data(D)))
    else:
        with open(fname, 'w') as outfile1:
            outfile1.write(json.dumps(json_graph.node_link_data(exo)))

'''
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
'''


def juliaJson2nx(tnet, dict, exogenous_G=False):
    d = {}
    for key in dict.keys():
        orig, dest = key.split(',')
        orig = orig.split("(")[1].replace('"', '').replace(' ', '')
        dest = dest.split(")")[0].replace('"', '').replace(' ', '')
        if "'" in orig:
            s = orig
        else:
            s = int(orig)
        if "'" in dest:
            t = dest
        else:
            t = int(dest)
        tnet.G_supergraph[s][t]['flow'] = dict[key]['flow']
        tnet.G_supergraph[s][t]['flowNoRebalancing'] = dict[key]['flow']
        if exogenous_G==False:
            tnet.G_supergraph[s][t]['t_k'] = travel_time(tnet,s,t, G_exo=exogenous_G)
        else:
            tnet.G_supergraph[s][t]['t_k'] = travel_time(tnet, s, t, G_exo=exogenous_G)


def solve_social_Julia(tnet, exogenous_G=False):
    # Save to json files
    nx2json(tnet.G, "tmp/G.json")
    nx2json(tnet.G_supergraph, "tmp/G_supergraph.json")
    nx2json(tnet.G_supergraph, "tmp/exogenous_G.json", exo=exogenous_G)

    js = json.dumps({str(k):v for k,v in tnet.g.items()})
    f = open("tmp/g.json", "w")
    f.write(js)
    f.close()

    f = open("tmp/fcoeffs.json", "w")
    f.write(str(tnet.fcoeffs))
    f.close()

    # Solve system-centric in julia
    shell("julia src/CARS.jl", printOut=True)

    # Parse results back
    dict_G = json2dict("tmp/out.json")
    juliaJson2nx(tnet, dict_G, exogenous_G=exogenous_G)
    # Get solve time
    f = open('tmp/solvetime.txt')
    line = f.readline()
    f.close()
    solvetime = float(line)
    shell("rm tmp/out.json", printOut=False)
    shell("rm tmp/G.json", printOut=False)
    shell("rm tmp/G_supergraph.json", printOut=False)
    shell("rm tmp/exogenous_G.json", printOut=False)
    shell("rm tmp/g.jsonn", printOut=False)
    shell("rm tmp/fcoeffs.json", printOut=False)
    shell("rm tmp/solvetime.txt", printOut=False)

    return solvetime
    #TODO: add delete funtion of out json




def solve_social_altruistic_Julia(tnet, exogenous_G=False):
    # Save to json files
    nx2json(tnet.G, "tmp/G.json")
    nx2json(tnet.G_supergraph, "tmp/G_supergraph.json")
    if exogenous_G != False:
        nx2json(exogenous_G, "tmp/exogenous_G.json", exo=True)
    else:
        nx2json(tnet.G, "tmp/exogenous_G.json", exo=False)

    js = json.dumps({str(k):v for k,v in tnet.g.items()})
    f = open("tmp/g.json", "w")
    f.write(js)
    f.close()

    f = open("tmp/fcoeffs.json", "w")
    f.write(str(tnet.fcoeffs))
    f.close()

    # Solve system-centric in julia
    shell("julia src/CARS_altruistic.jl", printOut=False)

    # Parse results back
    dict_G = json2dict("tmp/out.json")
    juliaJson2nx(tnet, dict_G)



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

def hist_flows(G, G_exo=True):
    if G_exo:
        norm_flows = [(G[i][j]['flow'] + G_exo[i][j]['flow']) / G[i][j]['capacity'] for i,j in G.edges()]
    else:
        norm_flows = [G[i][j]['flow'] / G[i][j]['capacity'] for i, j in G.edges()]
    #_ = plt.hist(norm_flows, bins='auto')
    #count, bins = np.histogram(norm_flows, bins=5)
    #print('bins:' + str(bins))
    print('max flow:' + str(round(max(norm_flows),2)))

    #fig, axs = plt.subplots(1, 1)
    #axs[0].hist(norm_flows, bins=5)
    #plt.show()
