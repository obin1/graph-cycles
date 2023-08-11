#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 13:59:32 2023

@author: psturm
"""

import numpy as np
import networkx as nx
from scipy.sparse import csr_matrix
from networkx.algorithms import bipartite
import random as rd




B = nx.read_edgelist("GCv14_Edgelist.csv",delimiter=", ")

nodelist = np.array(B.nodes)
nodebool = np.full(len(nodelist), True)

for r in range(1,914):
    Rstring = 'R'+str(r)
    Rloc = np.where(nodelist == Rstring)[0][0]
    nodebool[Rloc] = False
    

nodelist[nodebool]

    
G = bipartite.projected_graph(B, nodelist[nodebool], multigraph=True)
nx.write_edgelist(G, "GC_unipartite_edgelist.csv",  delimiter=',', data=False, encoding='utf-8')
