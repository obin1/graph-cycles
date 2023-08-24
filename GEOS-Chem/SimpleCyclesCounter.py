#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 24 11:37:46 2023

@author: psturm
"""


import numpy as np
import pandas as pd
# Note: tested with networkx 3.1, 2.8.4 did not have chordless cycles
import networkx as nx

# df = pd.read_csv("GCv14_Edgelist.csv")
df = pd.read_csv("/Users/psturm/Desktop/KPP-playground/GCv14/gckpp_EdgeList.csv")


# Create DiGraph
B = nx.DiGraph()

# Add edges. Note CSV file had an extra space after the comma 
B.add_edges_from([(df["from"][i],df["to"][i]) for i in range(0,len(df["to"]))])

# Create adjacency matrix 
Badj = nx.adjacency_matrix(B)
Badj = Badj.tocoo()
dt = {'row': Badj.row+1, 'col': Badj.col+1, 'data': Badj.data}
Badj2 = pd.DataFrame(data=dt)
Badj2.to_csv('GC14ReciprocityMatrixSparse.csv', index=False)

# Count simple, chordless cycles
reciprocal_cycles = list(nx.simple_cycles(B,length_bound=2))

# Note the length bound results in cycles of lower length as well
# Think of cycles_4 as "cycles up to length 4"
print("calculating 4-cycles")
cycles_4 = list(nx.simple_cycles(B,length_bound=4))
print("Removing 2 cycles from 4 cycles")
cycles_4 = [c for c in cycles_4 if len(c) == 4]

print("calculating 6-cycles")
cycles_6 = list(nx.simple_cycles(B,length_bound=6))
print("Removing 2,4 cycles from 6 cycles")
cycles_6 = [c for c in cycles_6 if len(c) == 6]

print("calculating 8-cycles")
cycles_8 = list(nx.simple_cycles(B,length_bound=8))
print("Removing 2,4,6 cycles from 8 cycles")
cycles_8 = [c for c in cycles_8 if len(c) == 8]



# # Now save these lists to csv files
print("Saving 4-cycles")
df_4 = pd.DataFrame(cycles_4)
df_4.to_csv("GCv14_4cycles_simple.csv",index=False)

print("Saving 6-cycles")
df_6 = pd.DataFrame(cycles_6)
df_6.to_csv("GCv14_6cycles_simple.csv",index=False)

print("Saving 8-cycles")
df_8 = pd.DataFrame(cycles_8)
df_8.to_csv("GCv14_8cycles_simple.csv",index=False)



