#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 10 16:32:49 2023

@author: psturm
"""

import numpy as np
import igraph as ig
from igraph import Graph
import matplotlib.pyplot as plt
import networkx as nx
from scipy.sparse import csr_matrix

# Link to small stratospheric mechanism: https://kpp.readthedocs.io/en/latest/getting_started/02_running_kpp_sample_mech.html
# Note bipartite tools are not imported into the networkx namespace
from networkx.algorithms import bipartite

# rows of sparse stoichiometric matrix from small_strato_StoichiomSP.f
rows = np.array([2,  2,  3,  2,  3,  2,  3,  1,  3,  1,  2,  
                 1,  3,  3,  4,  5,  2,  4,  5,  2,  4,  5 ])
# columns of sparse stoichiometric matrix from small_strato_StoichiomSP.f
cols = np.array([1,  2,  2,  3,  3,  4,  4,  5,  5,  6,  6,  
                 7,  7,  8,  8,  8,  9,  9,  9, 10, 10, 10 ])
# stoichiometric coefficients of sparse matrix from small_strato_StoichiomSP.f
vals =    np.array([ 2,   -1,    1,    1,   -1,   -1,   -1,    1,
                    -1,   -1,    1,   -1,   -1,   -1,   -1,    1,
                    -1,    1,   -1,    1,    1,   -1 ])


n_spc = np.max(rows)
n_rxn = np.max(cols)

species = ["O1D","O","O3","NO","NO2"]
reactions = ['R' + str(i) for i in range(1,np.max(cols)+1)]

elist = [(species[rows[i]-1], reactions[cols[i]-1], -vals[i]) for i in range(0,len(rows))]
for i in range(0,len(rows)):
    if elist[i][2]<0:
        elist[i] = (elist[i][1],elist[i][0],-elist[i][2])

# Generate weighted, directed incidence (or biadjacency) matrix 
# (also called the stoichiometric matrix)
# As in Sturm and Wexler (2022): https://doi.org/10.5194/gmd-15-3417-2022
# Note that networkx calls this a biadjacency matrix
A = np.zeros((5,10))
A[rows-1,cols-1] = vals
Asp = csr_matrix((vals, (rows-1, cols-1)))



# Note: 
# Create edge list of bipartite graph, with species first
# edges = [(rows[i]-1, n_spc+cols[i]-1) for i in range(0,len(rows))]
# g = ig.Graph.Bipartite( np.hstack(( [False]*n_spc, [True]*n_rxn)), edges,directed=True)
# g.es['weight'] = vals
# G = g.to_networkx()

# # 
# nx.incidence_matrix(G)

B = nx.DiGraph()
spc_nodes = np.unique(rows-1)
rxn_nodes = np.unique(cols)-1+n_spc

spc_list = ['O1D',',O,','O3','NO','NO2']
B.add_nodes_from(reactions, bipartite=0)
B.add_nodes_from(species, bipartite=1)
B.add_weighted_edges_from(elist)
color_map = ['salmon' if node in reactions else 'lightblue' for node in B]
size_map = [500 if node in reactions else 1000 for node in B]
size_map = [2000 if node == "O3" else 1000 for node in B]
nx.draw_kamada_kawai(B,node_color=color_map,node_size = size_map, with_labels=True)
plt.savefig('/Users/psturm/Desktop/Graph cycles/small_strato/small_strato.pdf', format='pdf')

# Project to unipartite species graph
# G = bipartite.projected_graph(B, spc_nodes,multigraph=True)

# Can we find  cycles?
# note there are other algorithms in this library
B_cycles = sorted(nx.simple_cycles(B))
# G_cycles = sorted(nx.simple_cycles(G))
