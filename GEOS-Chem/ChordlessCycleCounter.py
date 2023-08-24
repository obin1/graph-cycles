#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 10 22:55:34 2023

@author: psturm
"""

import numpy as np
import pandas as pd
# Note: tested with networkx 3.1, 2.8.4 did not have chordless cycles
import networkx as nx

df = pd.read_csv("GCv14_Edgelist.csv")
# df = pd.read_csv("/Users/psturm/Desktop/KGCv14_Edgelist.csv")


# Create DiGraph
B = nx.DiGraph()

# Add edges. Note CSV file had an extra space after the comma 
B.add_edges_from([(df["from"][i],df[" to"][i][1:]) for i in range(0,len(df[" to"]))])

# Count simple, chordless cycles
reciprocal_cycles = list(nx.chordless_cycles(B,length_bound=2))

# Note the length bound results in cycles of lower length as well
# Think of cycles_4 as "cycles up to length 4"
print("calculating 4-cycles")
cycles_4 = list(nx.chordless_cycles(B,length_bound=4))
print("calculating 6-cycles")
cycles_6 = list(nx.chordless_cycles(B,length_bound=6))
print("calculating 8-cycles")
cycles_8 = list(nx.chordless_cycles(B,length_bound=8))

# How many 8 cycles are there?
print("There are ", len(cycles_8)-len(cycles_6), " chordless 8 cycles")

# Remove the 2,4,6 cycles from the 8 cycles
print("Removing 2,4,6 cycles from 8 cycles")
cycles_8 = [c for c in cycles_8 if c not in cycles_6 and c not in cycles_4 and c not in reciprocal_cycles]

print("Removing 2,4 cycles from 6 cycles")
cycles_6 = [c for c in cycles_6 if c not in cycles_4 and c not in reciprocal_cycles]

print("Removing 2 cycles from 4 cycles")
cycles_4 = [c for c in cycles_4 if c not in reciprocal_cycles]

# Now save these lists to csv files
print("Saving 4-cycles")
df_4 = pd.DataFrame(cycles_4)
df_4.to_csv("GCv14_4cycles.csv",index=False)

print("Saving 6-cycles")
df_6 = pd.DataFrame(cycles_6)
df_6.to_csv("GCv14_6cycles.csv",index=False)

print("Saving 8-cycles")
df_8 = pd.DataFrame(cycles_8)
df_8.to_csv("GCv14_8cycles.csv",index=False)
