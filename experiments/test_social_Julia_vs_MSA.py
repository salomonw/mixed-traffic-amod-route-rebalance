import src.tnet as tnet
import src.CARS as cars
from src.utils import *
import numpy as np
import matplotlib.pyplot as plt
import copy

#netFile, gFile, fcoeffs, tstamp, dir_out = tnet.get_network_parameters('Braess1', experiment_name='Braess1_penRate_disjoint')
netFile, gFile, fcoeffs, tstamp, dir_out = tnet.get_network_parameters('EMA', experiment_name='EMA_system_centric_vs_altruistic')


g_multiplier_v = []
msa_v = []
total_msa_v=[]
nlp_v = []
total_nlp_v=[]
nlp_altruistic_v = []
total_nlp_altruisitc_v = []
print('------------------------------------------------------------------------------------')
print('Exog g \t Method \t AMoD tt \t Global tt')
print('------------------------------------------------------------------------------------')
for g_multiplier in np.linspace(0.025, 5, 20):

    for method in ['MSA', 'NLP', 'NLP_altruistic']: #['NLP', 'NLP_altruistic']:#[

        tNet = tnet.tNet(netFile=netFile, gFile=gFile, fcoeffs=fcoeffs)
        tNet.build_supergraph()
        pedestrian = [(u, v) for (u, v, d) in tNet.G_supergraph.edges(data=True) if d['type'] == 'p']
        connector = [(u, v) for (u, v, d) in tNet.G_supergraph.edges(data=True) if d['type'] == 'f']
        g_k = tnet.perturbDemandConstant(tNet.g, constant=g_multiplier)


        tNet_exog = copy.deepcopy(tNet)
        tNet_exog.set_g(g_k)
        tNet_exog.solveMSA()

        cars.G2supergraph(tNet_exog)


        if method == 'MSA':
            tNet.solveMSAsocial_supergraph(exogenous_G=tNet_exog.G_supergraph) #Revisar por que no deja evaluar el resultado
            msa = tnet.get_totalTravelTime(tNet.G_supergraph, fcoeffs, G_exogenous=tNet_exog.G_supergraph)
            tNet.G_supergraph = tNet_exog.G_supergraph.copy()
            cars.G2supergraph(tNet)
            msa_exo = tnet.get_totalTravelTime(tNet_exog.G_supergraph, fcoeffs, G_exogenous=tNet.G_supergraph)
            msa_v.append(msa/tNet.totalDemand)
            total_msa_v.append((msa+msa_exo)/(tNet.totalDemand+tNet_exog.totalDemand))
            print(str(g_multiplier) + '\t MSA \t' + str(msa/tNet.totalDemand) +'\t'+ str((msa+msa_exo)/(tNet.totalDemand+tNet_exog.totalDemand)))
        elif method == 'NLP':
            cars.solve_social_Julia(tNet, exogenous_G=tNet_exog.G_supergraph)
            nlp = tnet.get_totalTravelTime(tNet.G_supergraph, fcoeffs, G_exogenous=tNet_exog.G_supergraph)
            nlp_exo = tnet.get_totalTravelTime(tNet_exog.G_supergraph, fcoeffs, G_exogenous=tNet.G_supergraph)
            nlp_v.append(nlp/tNet.totalDemand)
            total_nlp_v.append((nlp+nlp_exo)/(tNet.totalDemand+tNet_exog.totalDemand))
            print(str(g_multiplier) + '\t NLP \t' + str(nlp/tNet.totalDemand) +'\t'+ str((nlp+nlp_exo)/(tNet.totalDemand+tNet_exog.totalDemand)))
        elif method == 'NLP_altruistic':
            cars.solve_social_altruistic_Julia(tNet, exogenous_G=tNet_exog.G_supergraph)
            nlp_altruistic = tnet.get_totalTravelTime(tNet.G_supergraph, fcoeffs, G_exogenous=tNet_exog.G_supergraph)
            nlp_altruistic_exo = tnet.get_totalTravelTime(tNet_exog.G_supergraph, fcoeffs, G_exogenous=tNet.G_supergraph)
            nlp_altruistic_v.append(nlp_altruistic/tNet.totalDemand)
            total_nlp_altruisitc_v.append((nlp_altruistic+nlp_altruistic_exo)/(tNet.totalDemand+tNet_exog.totalDemand))
            print(str(g_multiplier) + '\t NLPa \t' + str(nlp_altruistic/tNet.totalDemand) +'\t'+ str((nlp_altruistic+nlp_altruistic_exo)/(tNet.totalDemand+tNet_exog.totalDemand)))

    g_multiplier_v.append(g_multiplier)


fig, ax1 = plt.subplots()
ax1.plot(g_multiplier_v, nlp_v, '--',  label='NLP (AMoD)', color='orange')
ax1.plot(g_multiplier_v, msa_v, '--',  label='MSA (AMoD)', color='blue')
ax1.plot(g_multiplier_v, nlp_altruistic_v, '--',  label='NLP Altruistic (AMoD)', color='green')
ax1.plot(g_multiplier_v, total_nlp_v, '-',  label='NLP (Total)', color='orange')
ax1.plot(g_multiplier_v, total_msa_v, '-',  label='MSA (Total)', color='blue')
ax1.plot(g_multiplier_v, total_nlp_altruisitc_v, '-',  label='NLP Altruistic (Total)' , color='green')
plt.legend()
plt.xlabel('Exogenous vehicles demand multiplier')
plt.ylabel('Total travel time')

mkdir_n('results/' + dir_out)
plt.savefig('results/' + dir_out +'/costs.png', dpi=300)

from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
axins = zoomed_inset_axes(ax1, 30, loc=4) # zoom-factor: 2.5, location: upper-left
axins.plot(g_multiplier_v, nlp_v, '--',  label='NLP (AMoD)', color='orange')
axins.plot(g_multiplier_v, msa_v, '--',  label='MSA (AMoD)', color='blue')
axins.plot(g_multiplier_v, nlp_altruistic_v, '--',  label='NLP Altruistic (AMoD)', color='green')
axins.plot(g_multiplier_v, total_nlp_v, '-',  label='NLP (Total)', color='orange')
axins.plot(g_multiplier_v, total_msa_v, '-',  label='MSA (Total)', color='blue')
axins.plot(g_multiplier_v, total_nlp_altruisitc_v, '-',  label='NLP Altruistic (Total)' , color='green')
x1, x2, y1, y2 = 0.9 , 1.2 , 0.46 , 0.51 # specify the limits
axins.set_xlim(x1, x2) # apply the x-limits
axins.set_ylim(y1, y2) # apply the y-limits
plt.yticks(visible=False)
plt.xticks(visible=False)
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
mark_inset(ax1, axins, loc1=3, loc2=2, fc="none", ec="2.5")

mkdir_n('results/' + dir_out)
plt.savefig('results/' + dir_out +'/zoomed.png', dpi=300)

plt.show()