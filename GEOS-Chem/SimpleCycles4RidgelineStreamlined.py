#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 28 02:13:45 2023

@author: emywli
"""

import numpy as np
import networkx as nx
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})

# df = pd.read_csv("GCv14_Edgelist.csv")
df = pd.read_csv("gckpp_EdgeList.csv")
# Create DiGraph
B = nx.DiGraph()
# Add edges. Note CSV file had an extra space after the comma
B.add_edges_from([(df["from"][i],df["to"][i]) for i in range(0,len(df["to"]))])

# EdgeListFile = open('.csv','r')
# data = EdgeListFile.read()
# lines = data.split('\n')
# EdgeListFile.close()

# # eqnfile = open('/Users/emywli/Desktop/Research/Sources/EQN_NAMES Updated.txt', "rt")
# # data2 = eqnfile.read()
# # eqnlist = data2.split('\n')
# # eqnfile.close()

# obtaining species and converting to list from SPC_NAMES.txt
my_file = open('SPC_NAMES.txt', "rt")
data = my_file.read()
species = data.split()
my_file.close()

reactions = ['R' + str(i) for i in range(1,914)]

# # Create DiGraph
# B = nx.DiGraph()

# elistuw = []
# for line in range(0, len(lines)): 
#     temp = lines[line].split(', ')
#     elistuw.append(temp)

spc_list = species
# B.add_nodes_from(reactions, bipartite=0)
# B.add_nodes_from(species, bipartite=1)
# B.add_edges_from(elistuw)
# color_map = ['salmon' if node in reactions else 'lightblue' for node in B]
# size_map = [10 if node in reactions else 20 for node in B]

# Reading in stored cycles from csv 
cycles4 = open('GCv14_4cycles_simple.csv','r')
data4 = cycles4.read()
lines4 = data4.split('\n')
cycles4.close()

lines4.pop(0)
lines4.pop(-1)

cyclesfour = []
for i in range(0, len(lines4)): 
    temp = lines4[i].split(',')
    cyclesfour.append(temp)
    
totalcycles = cyclesfour

# list to store the distribution data 
timescalelin = []

# list to store the repeat location names corresponding to the number of timescale values for that location
tiles = []

# location names as they are in the file names, in order of how we want them to show up on graph
locations = [['AtlanticOcean', '1500'], ['IndianOcean', '0600'], ['PacificOcean', '1900'], ['Graciosa', '1200'], ['CapeGrim', '0200'], ['ElDjouf', '1200'], ['Amazon', '1600'], ['Borneo', '0400'], ['Congo', '1100'], ['Ozarks', '1700'], ['Beijing', '0400'], ['Kinshasa', '1100'], ['LosAngeles', '1900'], ['Paris', '1000'], ['Utqiagvik', '2000'], ['McMurdo', '0000']]
# locations = [['Amazon', '1600']]
# location names as how we want them to appear on the graph
graphlocations = ['Atlantic Ocean', 'Indian Ocean', 'Pacific Ocean', 'Graciosa', 'Kennaook', 'El Djouf', 'Amazon', 'Borneo', 'Congo', 'Ozarks', 'Beijing', 'Kinshasa', 'Los Angeles', 'Paris', 'Utqiagvik', 'McMurdo Station']
# filename = 'L1 Noon Samples/AtlanticOcean_L1_20180702_2000.txt'
for location in locations: 
    filename = 'L1 Noon Samples/' + str(location[0]) + '_L1_20180702_' + str(location[1]) + '.txt'
    sample_file = open(filename, "rt")
    sampledata = sample_file.read()
    sample = sampledata.split('\n ')
    sample_file.close()
    # sampleconcentrations = []
    # for c in range(0, len(sample)): 
    #     if sample[c][0] == 'C': 
    #         temp = sample[c].split()
    #         sampleconcentrations.append(temp[4])
    # samplerateconst = []
    # for k in range(0, len(sample)): 
    #     if sample[k][0] == 'R': 
    #         temp = sample[k].split()
    #         samplerateconst.append(temp[4])
    rates = [] 
    for r in range(0, len(sample)): 
        if sample[r][0] == 'A': 
            temp = sample[r].split()
            rates.append(temp[4])    
    for i in range(0, len(rates)):
        rates[i] = float(rates[i])
    # list to hold residence times for each cycle in totalcycles
    cyclest = []
    # # track = 0
    for i in range(0, len(totalcycles)): 
        rnums = []
        sub_cyclest = 0
        # if the cycle on current iteration starts with a reaction or species
        if totalcycles[i][0] in reactions:
            for j in range(0, len(totalcycles[i]), 2):
                rnums.append(totalcycles[i][j].replace('R', ''))
        else: 
            for j in range(1, len(totalcycles[i]), 2):
                rnums.append(totalcycles[i][j].replace('R', ''))
        # add time per reaction to get residence time for the full cycle on current iteration 
        for m in range(0, len(rnums)): 
            if rates[int(rnums[m])-1] == float(0):
                sub_cyclest += float('inf')
            else: 
                sub_cyclest += 1/rates[int(rnums[m])-1]
        cyclest.append([sub_cyclest, totalcycles[i]])
    # # to save a copy of the residence times list before using the original cyclest one in sorting cycles
    # cycleststore = []
    # for i in range(0, len(cyclest)):
    #     cycleststore.append(cyclest[i])
    # # for the histogram
    # # removing the residence times that are 'inf' in order to be able to calculate bin size, noting that inf would skew it higher
    cyclest_noninf = []
    for i in range(0, len(cyclest)):
        if cyclest[i][0] != float('inf'): 
            cyclest_noninf.append([cyclest[i][0], cyclest[i][1]])
    cyclest_noninf.sort()
    # do we want to store the version of cyclest_noninf with both the timescale and their corresponding cycle? / helps with analysis? 
    times = []
    for i in range(0, len(cyclest_noninf)): 
        times.append(float(cyclest_noninf[i][0]))
    timescalelin += times
    # getting index of current iteration's location to put into tiles list
    index = locations.index(location)
    tiles += len(times) * [graphlocations[index]]
    # df_4 = pd.DataFrame(cyclest_noninf)
    # df_4.to_csv('/Users/emywli/Desktop/Research/L56 Noon Cycles/AtlanticOcean_L56_4cyclestimes_simple.csv',index=False)
    
timescalelin = np.array(timescalelin)
timescale = np.log10(timescalelin)
df = pd.DataFrame(dict(timescale=timescale, tiles=tiles))
# m = df.tiles.map(ord)
# df["timescale"] += m

# Initialize the FacetGrid object
# pal = sns.cubehelix_palette(3, rot=-.25, light=.7)
# pal = sns.color_palette(palette='magma_r', n_colors=3)
pal = sns.cubehelix_palette(16, rot=-.25, light=.7)
pal[15] = pal[0]
pal[14] = pal[1]
pal[0] = pal[1] = pal[2] = pal[3] = pal[4] = [0.2154516905067, 0.295, 0.5078367867510504]
pal[5] = [0.78, 0.51, 0.1892045842]
pal[6] = pal[7] = pal[8] = pal[9] = [0.12071162840208301, 0.4357915995132193, 0.2463679091477368]
pal[10] = pal[11] = pal[12] = pal[13] = [0.7, 0.7, 0.65]

# pal[13] = [0.78, 0.51, 0.1892045842]
tiles = sns.FacetGrid(df, row="tiles", hue="tiles", aspect=15, height=.5, palette=pal)

# Draw the densities in a few steps
tiles.map(sns.kdeplot, "timescale",
      bw_adjust=.5, clip_on=False,
      fill=True, alpha=1, linewidth=1.5)
tiles.map(sns.kdeplot, "timescale", clip_on=False, color="w", lw=2, bw_adjust=.5)

# passing color=None to refline() uses the hue mapping
tiles.refline(y=0, linewidth=2, linestyle="-", color=None, clip_on=False)


# Define and use a simple function to label the plot in axes coordinates
def label(x, color, label):
    ax = plt.gca()
    ax.text(0, .2, label, fontweight="bold", color=color,
            ha="left", va="center", transform=ax.transAxes)
    # ax.canvas.draw()
    # plt.tight_layout()

tiles.map(label, "timescale")

# Set the subplots to overlap
tiles.figure.subplots_adjust(hspace=-0.15)

# Remove axes details that don't play well with overlap
tiles.set_titles("")
tiles.set(yticks=[], ylabel="")
tiles.despine(bottom=True, left=True)

# plt.savefig('Noon Plots/Timescale Distributions/RidgelineLocationsL1NoonGrouped.pdf', bbox_inches="tight")
