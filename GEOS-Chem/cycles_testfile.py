#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 10:22:55 2023

@author: emywli
"""

import numpy as np
import networkx as nx
from scipy.sparse import csr_matrix
from networkx.algorithms import bipartite
import random as rd

rd.seed(42)

# Pulling biadjaceny matrix info from the A.csv file
Acsvfile = open('/Users/emywli/Desktop/Research/Sources/A.csv','r')
data = Acsvfile.read()
Adata = data.split()
# remove the extra commas between value entries in the A.csv file
Adatalist = []
for i in range(0, len(Adata)):
    if Adata[i] != ',':
        Adatalist.append(float(Adata[i])) # how to get more precision in the floating value?
Acsvfile.close()

# for i in range(0, len(Adatalist)-1): 
#     int(Adatalist[i])


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


n_spc = max(rows) # n_spc returning 287 when there are 291 species in txt file? 
n_rxn = max(cols)

# obtaining species and converting to list from SPC_NAMES.txt
my_file = open('/Users/emywli/Desktop/Research/Sources/SPC_NAMES.txt', "rt")
data = my_file.read()
species = data.split()
my_file.close()
reactions = ['R' + str(i) for i in range(1,int(n_rxn)+1)]

elist = [(species[int(rows[i])-1], reactions[int(cols[i])-1], -vals[i]) for i in range(0,len(rows))]
for i in range(0,len(rows)):
    if elist[i][2]<0:
        elist[i] = (elist[i][1],elist[i][0],-elist[i][2])

# Generate weighted, directed incidence (or biadjacency) matrix 
# (also called the stoichiometric matrix)
# As in Sturm and Wexler (2021): https://doi.org/10.5194/gmd-15-3417-2022
# Note that networkx calls this a biadjacency matrix

rowsints = [int(num) for num in rows]
colsints = [int(num) for num in cols]
rowsintsarr = np.array(rowsints)
colsintsarr = np.array(colsints)

A = np.zeros((int(n_spc), int(n_rxn)))
A[rowsintsarr-1,colsintsarr-1] = vals
Asp = csr_matrix((vals, (rowsintsarr-1, colsintsarr-1)))

# seed
# for now, randomized list of reaction rates (to serve as the list "r") for all 913 reactions
randomrates = []
for i in range(0, int(n_rxn)):
    ratei = rd.random()
    randomrates.append(ratei)
    

B = nx.DiGraph()
spc_nodes = np.unique(rowsintsarr-1)
rxn_nodes = np.unique(colsintsarr)-1+int(n_spc)

spc_list = species # currently not being used but could be used in creation of G (instead of "species" parameter)
B.add_nodes_from(reactions, bipartite=0)
B.add_nodes_from(species, bipartite=1)
B.add_weighted_edges_from(elist)
# color_map = ['salmon' if node in reactions else 'lightblue' for node in B]
# size_map = [10 if node in reactions else 20 for node in B]
# kamada takes ~a million years for this 
# nx.draw_kamada_kawai(B,node_color=color_map,node_size = size_map, width=0.2, with_labels=True)
# nx.draw(B,node_color=color_map,node_size = size_map, width = 0.2, with_labels=True)

# Project to unipartite species graph
# G = bipartite.projected_graph(B, species, multigraph=True)
# nx.draw(G, node_size = 10, width = 0.2, with_labels=True)

edgesrB = list(B.edges)
# splitting the R[num] string i.e. G.edges[i][2] since reaction numbers are no longer 1 digit -- need to remove R from the string
rxnnums = []
for edge in B.edges:
    x = (edge[1]).replace('R', '')
    rxnnums.append(x)
    
# t_edges = np.array([ randomrates[int(rxnnums[i])-1] for i in range(0, len(rxnnums))])
# simple_edges = list(nx.Graph(G).edges())
# t_simple = np.zeros(len(simple_edges))
# for i in range(0,len(t_simple)):
#     for edge in G.edges:
#         if set(edge[0:2]) == set(simple_edges[i]):
#             if edge[0:2] == simple_edges[i]:
#                 r_simple[i] += 1/randomrates[int(edge[2][1])-1]  
            
cycleslist = []    
# seems like this generator route still makes the most sense instead of looping over all 4468 rows of A and finding all the cycles that start with each element
B_cycles = nx.simple_cycles(B)
# for i in range(0,1000):
#     cycleslist.append(next(B_cycles))  

num_cycle = 0
x = []
while num_cycle <= 100000000:
    # cycleslist.append(next(B_cycles))
    x = next(B_cycles)
    num_cycle += 1
# for num in B_cycles: 
#     cycleslist.append(num)

# cyclesarray = np.array(cycleslist, dtype=object)
# cyclesset = set(cycleslist)
# cyclesunique = list(cyclesset)

# for i in range(0, len(cycleslist)): 
#     if cycleslist[i] not in cyclesunique: 
#         cyclesunique.append(cycleslist[i])
        
    



