#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 14:07:46 2023

@author: emywli
"""

import numpy as np
import networkx as nx
from scipy.sparse import csr_matrix
from networkx.algorithms import bipartite
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap, Normalize
import random as rd

# function to pull subgraph, taking in a list of selected species, from GC14 and plot it 

Acsvfile = open('A.csv','r')
data = Acsvfile.read()
Adata = data.split()
Adatalist = []
for i in range(0, len(Adata)):
    if Adata[i] != ',':
        Adatalist.append(float(Adata[i])) # how to get more precision in the floating value?
Acsvfile.close()


# rows of sparse stoichiometric matrix 
rowslist = []
for i in range(0, len(Adatalist)-1, 3): 
    rowslist.append(Adatalist[i])
    
rows = np.array(rowslist)

# columns of sparse stoichiometric matrix 
colslist = []
for i in range(1, len(Adatalist), 3): 
    colslist.append(Adatalist[i])
    
cols = np.array(colslist)

# stoichiometric coefficients of sparse matrix 
valslist = []
for i in range(2, len(Adatalist), 3): 
    valslist.append(Adatalist[i])
    
vals = np.array(valslist)


n_spc = max(rows) 
n_rxn = max(cols)

# obtaining species and converting to list from SPC_NAMES.txt
my_file = open('SPC_NAMES.txt', "rt")
data = my_file.read()
species = data.split()
my_file.close()
reactions = ['R' + str(i) for i in range(1,int(n_rxn)+1)]

elist = [(species[int(rows[i])-1], reactions[int(cols[i])-1], -vals[i]) for i in range(0,len(rows))]
for i in range(0,len(rows)):
    if elist[i][2]<0:
        elist[i] = (elist[i][1],elist[i][0],-elist[i][2])

rowsints = [int(num) for num in rows]
colsints = [int(num) for num in cols]
rowsintsarr = np.array(rowsints)
colsintsarr = np.array(colsints)

A = np.zeros((int(n_spc), int(n_rxn)))
A[rowsintsarr-1,colsintsarr-1] = vals
Asp = csr_matrix((vals, (rowsintsarr-1, colsintsarr-1)))


B = nx.DiGraph()
spc_nodes = np.unique(rowsintsarr-1)
rxn_nodes = np.unique(colsintsarr)-1+int(n_spc)

spc_list = species
B.add_nodes_from(reactions, bipartite=0)
B.add_nodes_from(species, bipartite=1)
B.add_weighted_edges_from(elist)
color_map = ['salmon' if node in reactions else 'lightblue' for node in B]
size_map = [10 if node in reactions else 20 for node in B]
# nx.draw(B,node_color=color_map,node_size = size_map, width = 0.2, with_labels=True)


def subgraph(selectedspclist): 
    edgelist = []
    rxnnodes = []
    for spc in selectedspclist: 
        spcindex = species.index(spc)
        for i in range(0, len(rows)): 
            if int(rows[i]) == spcindex+1: 
                edgelist.append((spc, reactions[int(cols[i])-1], -vals[i])) 
                # adding to list the reaction numbers that each species is involved in
                rxnnodes.append(reactions[int(cols[i])-1])
    # print(rxnnodes)
    for i in range(0,len(edgelist)):
        if edgelist[i][2]<0:
            edgelist[i] = (edgelist[i][1],edgelist[i][0],-edgelist[i][2])
    S = nx.DiGraph()
    S.add_nodes_from(selectedspclist, bipartite=1)
    S.add_nodes_from(rxnnodes, bipartite = 0)
    S.add_weighted_edges_from(edgelist)
    color_map = ['salmon' if node in reactions else 'lightblue' for node in S]
    size_map = [100 if node in reactions else 200 for node in S]
    # nx.draw_kamada_kawai(S,node_color=color_map,node_size = size_map, width = 0.4, with_labels=True)
    nx.draw(S,node_color=color_map,node_size = size_map, width = 0.4, with_labels=True)
    return edgelist       
    

# this version of the function takes an additional input, the minimum number of species out of the specified list that are required to be present
# in a reaction for it to be included/plotted in the subgraph 
def subgraph2(selectedspclist, min_spccount):
    edgelist = []
    rxnnodes = []
    if min_spccount > len(selectedspclist):
        print("Note: Need minimum required-number-of-present-species input to be less than or equal to the size of list of species given.")
    for spc in selectedspclist: 
        spcindex = species.index(spc)
        for i in range(0, len(rows)): 
            if int(rows[i]) == spcindex+1: 
                # get the reaction number the species is involved in
                Rnum = cols[i]
                # to track how many species in the selectedspclist input are in this reaction
                spccount = 1
                # # check for the indices <= 20 first because the method under following "else:" block for checking 20 indices before and after won't work if current index is <20
                # # more robust way than running through all 4468 each time but for now will just be doing that
                # if i <= 20: 
                for j in range(0, len(rows)):
                    if cols[j] == Rnum and j != i: 
                        if species[int(rows[j])-1] in selectedspclist:
                            spccount += 1   
                # else: 
                #     # check 20 rows before and after current reaction to get all species involved in the reaction 
                #     # and check whether it is in the input selected species list (though 17 before/after should suffice since mode is 17, for rxn 536 with most species involved)
                #     for j in range(i-20, i+20): 
                #         if cols[j] == Rnum: 
                #             if species[int(rows[j])-1] in selectedspclist:
                #                 spccount += 1  
                # print(spccount, spc)
                if spccount >= min_spccount:
                    # if int(cols[i]) != 231 and int(cols[i]) != 224:
                    edgelist.append((spc, reactions[int(cols[i])-1], -vals[i])) 
                    # adding to list the reaction numbers that each species is involved in to make node list for subgraph
                    rxnnodes.append(reactions[int(cols[i])-1])
                    # rxnnodes.append(reactions[262])
                    # print(rxnnodes)
    # going back through and adjusting direction of arrow based on whether each edge's species component is a reactant or product
    for i in range(0,len(edgelist)):
        if edgelist[i][2]<0:
            edgelist[i] = (edgelist[i][1],edgelist[i][0],-edgelist[i][2])
    S = nx.DiGraph()
    S.add_nodes_from(selectedspclist, bipartite=1)
    S.add_nodes_from(rxnnodes, bipartite = 0)
    S.add_weighted_edges_from(edgelist)
    color_map = ['salmon' if node in reactions else 'lightblue' for node in S]
    size_map = [200 if node in reactions else 400 for node in S]
    # nx.draw_kamada_kawai(S,node_color=color_map,node_size = size_map, width = 0.4, with_labels=True)
    nx.draw(S,node_color=color_map,node_size = size_map, width = 0.5, with_labels=True)
    return edgelist       


# general subgraph based on selected species and reactions desired
listspc = ['ClO', 'Cl', 'O2', 'O','O3', 'HCl', 'MO2', 'H2O', 'OH', 'CH4']
listrxn = ['R249', 'R261', 'R263', 'R268']
H = B.subgraph(listspc+listrxn)

# subgraph by element
listelem = ['I']
listelemrxn = []
elemnodes = []
elemnodesrxn = set([])
for i in range(0, len(listelem)): 
    for j in range(0, len(species)): 
        if listelem[i] in species[j]:
            elemnodes.append(species[j])
            for spcnum in range(0, len(rows)): 
                if int(rows[spcnum])-1 == j:
                    elemnodesrxn.add(reactions[int(cols[spcnum])-1])
# H = B.subgraph(elemnodes+list(elemnodesrxn))
# nx.draw_kamada_kawai(H, with_labels=True)


# list of rates for reactions 249, 261, 263, 268 (from sample IndianOcean_L10_20190702_1200.txt)
# rateslist = [('R249', 3690.71620231313), ('R261', 10126.0720521017), ('R263', 913911.1158992073), ('R268', 0.290583717584288)]
# rateslist = [('R249', 1.475685092915952e-97), ('R261', 2.180329006503819e-97), ('R263', 1.519736987858295e-97), ('R268',  1.229884761588709e-97)]

# new rates from McMurdo_L35_20180930_2100.txt from Obin
# rateslist = [('R249', 711.349599979795), ('R261', 3567.26843635985), ('R263', 942740.565436074), ('R268', 3253.40633883372)]
rateslist = [('R249', 1), ('R261', 200000), ('R263', 500000), ('R268', 800000)]

# for i in range(0, rxnnums): 
rxnrates = []
for edge in H.edges: 
    for i in range(0, len(rateslist)): 
        if edge[0] == rateslist[i][0] or edge[1] == rateslist[i][0]: 
            rxnrates.append(rateslist[i][1])
            # rxnrates.append(rateslist[i][1])

rxnrates = np.array(rxnrates)
node_sizes = 30*np.log10(1e9) # just so every decade is an order of magnitude
M = H.number_of_edges()
edge_colors = rxnrates # r/(np.max(r)+1e-2)# range(2, 35)# 
# edge_colors = [10**11, 20, 10**9, 2000, 2000000, 2000000000, 1000, 10000, 100000, 10000000, 30, 50000000000000000]
# to get residence times as edge colors rather than pure rates
# edge_colors = np.divide(1, rxnrates) 
# edge_alphas =  1+np.log10(rxnrates)/5
# edge_alphas =  np.divide(1, rxnrates)
# edge_alphas =  1+(rxnrates)/5
cmap = plt.cm.plasma  



seed = 93
pos = nx.random_layout(H, seed=seed)
pos["R249"] += (-0.4, -0.1)
pos["R261"] += (0.7, -0.2)
pos["R263"] += (0.7, -0.1)
pos["R268"] += (0.5, -0.1)
pos["ClO"] += (0.1, 0.05)
pos["CH4"] += (-0.4, -0.05)
pos["Cl"] += (0.1, 0)
pos["O2"] += (0.7, -0.3)
pos["O3"] += (0.5, 0.49)
pos["O"] += (-0.1, 0.00)
# pos["HCl"] += (-0.15, 0.2)
pos["MO2"] += (0.15, 0.2)
# pos["H2O"] += (-0.4, -0.2)
pos["OH"] += (-0.48, 0.4)


nodes = nx.draw_networkx_nodes(H, pos, node_size=node_sizes, node_color="lightblue")
nx.draw_networkx_labels(H, pos,font_size=6)

edge_color = rxnrates
norm = mpl.colors.Normalize(vmin=np.min(edge_colors), vmax=np.max(edge_colors))
pc = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
pc.set_array(edge_colors)

# nx.draw(MDG, pos, with_labels=True, edge_color = (1,1,1))
# for e in G.edges:
# edgesdrawn = []
for e, ec in zip(H.edges, edge_color):
    color = cmap(norm(ec))
    plt.gca().annotate("",
                       xy=pos[e[1]], 
                       xycoords='data',
                       xytext=pos[e[0]], 
                       textcoords='data',
                       arrowprops=dict(arrowstyle="->", color=color, # edge_color=edge_colors, 
                                        # edge_cmap=cmap,
                                       shrinkA=10, shrinkB=10,
                                       patchA=None, patchB=None,
                                       connectionstyle="arc3,rad=rrr".replace('rrr',str(rd.random()*0.5+0.1)))
                      )
    # edgesdrawn.append(e)


# pc = plt.cm.ScalarMappable(cmap=cmap)
# pc.set_clim(vmin=np.min(edge_colors), vmax=np.max(edge_colors))
# pc = mpl.collections.PatchCollection(e, cmap, norm=mpl.colors.LogNorm())
# pc = plt.cm.ScalarMappable(cmap=cmap, norm=mpl.colors.LogNorm())


# plt.xlabel('text')
# plt.gca().xaxis.xlabel('time (s)')
# plt.gca().xaxis.set_units('inch')
# pc.set_label('seconds/(molec/cm^3)')

plt.axis('off')
cbar = plt.colorbar(pc)
cbar.set_label('seconds/(molec/cm^3)', rotation=270, labelpad=18) 
plt.colorbar(pc)
# set_label('seconds/(molec/cm^3)', rotation=270, labelpad=18)
# plt.tight_layout()
# plt.show()
# plt.savefig('/Users/emywli/Desktop/Research/MolinaRowlandViz9.pdf')
    
    
    
    
    
    
    
