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
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})

# Ridgeline plotting reference: https://seaborn.pydata.org/examples/kde_ridgeplot.html

df = pd.read_csv("gckpp_EdgeList.csv")
# Create DiGraph
B = nx.DiGraph()
# Add edges. Note CSV file had an extra space after the comma
B.add_edges_from([(df["from"][i],df["to"][i]) for i in range(0,len(df["to"]))])

# obtaining species and converting to list from SPC_NAMES.txt
my_file = open('SPC_NAMES.txt', "rt")
data = my_file.read()
species = data.split()
my_file.close()

reactions = ['R' + str(i) for i in range(1,914)]

# Reading in stored cycles from csv 
cycles = open('GCv14_4cycles_simple.csv','r')
data = cycles.read()
lines = data.split('\n')
cycles.close()

lines.pop(0)
lines.pop(-1)

cycleslist = []
for i in range(0, len(lines)): 
    temp = lines[i].split(',')
    cycleslist.append(temp)
    
totalcycles = cycleslist

# list to store the data for distributions
timescalelin = []

# list to store the repeat location names corresponding to the number of timescale values for that location
tiles = []

# list to store the longest timescales + corresponding cycle in each location 
fastestcycles = []
# location names as they are in the file names, in order of how we want them to show up on graph
locations = [['AtlanticOcean', '1500'], ['IndianOcean', '0600'], ['PacificOcean', '1900'], ['Graciosa', '1200'], ['CapeGrim', '0200'], ['ElDjouf', '1200'], ['Amazon', '1600'], ['Borneo', '0400'], ['Congo', '1100'], ['Ozarks', '1700'], ['Beijing', '0400'], ['Kinshasa', '1100'], ['LosAngeles', '1900'], ['Paris', '1000'], ['Utqiagvik', '2000'], ['McMurdo', '0000']]
# location names as how we want them to appear on the graph
graphlocations = ['Atlantic Ocean', 'Indian Ocean', 'Pacific Ocean', 'Graciosa', 'Kennaook', 'El Djouf', 'Amazon', 'Borneo', 'Congo', 'Ozarks', 'Beijing', 'Kinshasa', 'Los Angeles', 'Paris', 'Utqiagvik', 'McMurdo Station']
for location in locations: 
    filename = 'L1 Noon Samples/' + str(location[0]) + '_L1_20180702_' + str(location[1]) + '.txt'
    sample_file = open(filename, "rt")
    sampledata = sample_file.read()
    sample = sampledata.split('\n ')
    sample_file.close()
    rates = [] 
    for r in range(0, len(sample)): 
        if sample[r][0] == 'A': 
            temp = sample[r].split()
            rates.append(temp[4])    
    for i in range(0, len(rates)):
        rates[i] = float(rates[i])
    # list to hold timescales for every cycle in totalcycles
    cyclest = []
    # # track = 0
    for i in range(0, len(totalcycles)): 
        rnums = []
        # sub_cyclest will contain, in order: full cycle timescale, timescale of longest reaction, ratio of longest reaction to full cycle timescale, the spcs/rxns of cycle
        sub_cyclest = 0
        # track longest timescale after dividing 1/(slowest rate) below 
        max_timescale = 0.0
        # if the cycle on current iteration starts with a reaction or species
        if totalcycles[i][0] in reactions:
            for j in range(0, len(totalcycles[i]), 2):
                rnums.append(totalcycles[i][j].replace('R', ''))
        else: 
            for j in range(1, len(totalcycles[i]), 2):
                rnums.append(totalcycles[i][j].replace('R', ''))
        # add time per reaction to get full cycle timescale on current iteration 
        for m in range(0, len(rnums)): 
            if rates[int(rnums[m])-1] == float(0):
                sub_cyclest += float('inf')
            else: 
                val = 1/rates[int(rnums[m])-1]
                sub_cyclest += val
                if val > max_timescale: 
                    max_timescale = val
        longest_to_full = max_timescale / sub_cyclest
        cyclest.append([sub_cyclest, max_timescale, longest_to_full, totalcycles[i]])
    # removing the cycle timescales that are 'inf' for realistic analysis, noting that inf would skew distributions higher
    cyclest_noninf = []
    for i in range(0, len(cyclest)):
        if np.log10(cyclest[i][0]) <= 20:
            cyclest_noninf.append([cyclest[i][0], cyclest[i][1], cyclest[i][2], cyclest[i][3]])
    cyclest_noninf.sort()
    fastestcycles.append([location[0], cyclest_noninf[0][3], cyclest_noninf[0][0], cyclest_noninf[0][1], cyclest_noninf[0][2]])
    times = []
    for i in range(0, len(cyclest_noninf)): 
        times.append(float(cyclest_noninf[i][0]))
    timescalelin += times
    # getting index of current iteration's location to put into tiles list
    index = locations.index(location)
    tiles += len(times) * [graphlocations[index]]
    
df_fast = pd.DataFrame(fastestcycles)
# df_fast.to_csv('L1NoonFastestCycles.csv',index=False)

timescalelin = np.array(timescalelin)
timescale = np.log10(timescalelin)
df = pd.DataFrame(dict(timescale=timescale, tiles=tiles))

# color palette, each category of location having a different color
pal = sns.cubehelix_palette(16, rot=-.25, light=.7)
pal[15] = pal[0]
pal[14] = pal[1]
pal[0] = pal[1] = pal[2] = pal[3] = pal[4] = [0.2154516905067, 0.295, 0.5078367867510504]
pal[5] = [0.78, 0.51, 0.1892045842]
pal[6] = pal[7] = pal[8] = pal[9] = [0.12071162840208301, 0.4357915995132193, 0.2463679091477368]
pal[10] = pal[11] = pal[12] = pal[13] = [0.7, 0.7, 0.65]

tiles = sns.FacetGrid(df, row="tiles", hue="tiles", aspect=15, height=.5, palette=pal)

# mapping the densities
tiles.map(sns.kdeplot, "timescale",
      bw_adjust=.5, clip_on=False,
      fill=True, alpha=1, linewidth=1.5)
tiles.map(sns.kdeplot, "timescale", clip_on=False, color="w", lw=2, bw_adjust=.5)

# passing color=None to refline() uses the hue mapping
tiles.refline(y=0, linewidth=2, linestyle="-", color=None, clip_on=False)

def label(x, color, label):
    ax = plt.gca()
    ax.text(0, .2, label, fontweight="bold", color=color,
            ha="left", va="center", transform=ax.transAxes)
    ax.set_xlim(-12, 22)


tiles.map(label, "timescale")

# Set the subplots to overlap
tiles.figure.subplots_adjust(hspace=-0.15)

    
# Remove axes details that don't play well with overlap
tiles.set_titles("")

ticks = plt.xticks()
xtick_string = [rf'$10^{str({i})}$' for i in range(-10, 25, 5)]

tiles.set(xticks=range(-10, 25, 5))
tiles.set_xticklabels(xtick_string)


tiles.set(xlabel='Timescale (seconds/cycle)')
tiles.set(yticks=[], ylabel="")

tiles.despine(bottom=True, left=True)

plt.show
# plt.savefig('/Users/emywli/Desktop/Research/RidgelineL1Noon_Manuscript.png', dpi=300, bbox_inches="tight")

