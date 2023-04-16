#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 22:26:50 2023

@author: emywli
"""

import numpy as np
import networkx as nx
from networkx.algorithms import bipartite


# Pulling biadjaceny matrix info from the A.csv file
Acsvfile = open('BiadjacencyMatrix.csv','r')
data = Acsvfile.read()
Adata = data.split()
# remove the extra commas between value entries in the A.csv file
Adatalist = []
for i in range(6, len(Adata)):
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

n_spc = max(rows) # n_spc returning 287 when there are 291 species in txt file? 
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
    

B = nx.DiGraph()
spc_nodes = np.unique(rowsintsarr-1)
rxn_nodes = np.unique(colsintsarr)-1+int(n_spc)

B.add_nodes_from(reactions, bipartite=0)
B.add_nodes_from(species, bipartite=1)
B.add_weighted_edges_from(elist)

cycleslist = []    
B_cycles = nx.simple_cycles(B)

num_cycle = 0
while num_cycle < 1e9:
    next(B_cycles)
    num_cycle += 1
