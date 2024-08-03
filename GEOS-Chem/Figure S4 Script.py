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

# obtaining species and converting to list from SPC_NAMES.txt
my_file = open('SPC_NAMES.txt', "rt")
data = my_file.read()
species = data.split()
my_file.close()

reactions = ['R' + str(i) for i in range(1,914)]

# Reading in stored cycles from csv 
cyclesstored = open('GCv14_4cycles_simple.csv','r')
data = cyclesstored.read()
lines = data.split('\n')
cyclesstored.close()

lines.pop(0)
lines.pop(-1)

cycleslist = []
for i in range(0, len(lines)): 
    temp = lines[i].split(',')
    cycleslist.append(temp)
    
totalcycles = cycleslist
OHcycles = []
nonOHcycles = []
for i in range(0, len(totalcycles)):
    for j in range(0, len(totalcycles[i])): 
        if totalcycles[i][j] == 'OH': 
            OHcycles.append(totalcycles[i])
            break
        elif j == len(totalcycles[i])-1:
            nonOHcycles.append(totalcycles[i])

            
# list to store the distribution data 
timescalelin = []

# list to store the repeat location names corresponding to the number of timescale values for that location
tiles = []

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
    # list to hold residence times for each cycle in OHcycles or nonOHcycles
    # can replace OHcycles with nonOHcycles and vice versa below, to plot the cycles involving OH and not involving OH, respectively
    cyclest = []
    for i in range(0, len(nonOHcycles)): 
        rnums = []
        sub_cyclest = 0
        # if the cycle on current iteration starts with a reaction or species
        if nonOHcycles[i][0] in reactions:
            for j in range(0, len(nonOHcycles[i]), 2):
                rnums.append(nonOHcycles[i][j].replace('R', ''))
        else: 
            for j in range(1, len(nonOHcycles[i]), 2):
                rnums.append(nonOHcycles[i][j].replace('R', ''))
        # add time per reaction to get residence time for the full cycle on current iteration 
        for m in range(0, len(rnums)): 
            if rates[int(rnums[m])-1] == float(0):
                sub_cyclest += float('inf')
            else: 
                sub_cyclest += 1/rates[int(rnums[m])-1]
        cyclest.append([sub_cyclest, nonOHcycles[i]])
    cyclest_noninf = []
    for i in range(0, len(cyclest)):
        if np.log10(cyclest[i][0]) <= 20:
            cyclest_noninf.append([cyclest[i][0], cyclest[i][1]])
    cyclest_noninf.sort()
    times = []
    for i in range(0, len(cyclest_noninf)): 
        times.append(float(cyclest_noninf[i][0]))
    timescalelin += times
    # getting index of current iteration's location to put into tiles list
    index = locations.index(location)
    tiles += len(times) * [graphlocations[index]]
  
timescalelin = np.array(timescalelin)
timescale = np.log10(timescalelin)
df = pd.DataFrame(dict(timescale=timescale, tiles=tiles))

pal = sns.cubehelix_palette(16, rot=-.25, light=.7)
pal[15] = pal[0]
pal[14] = pal[1]
pal[0] = pal[1] = pal[2] = pal[3] = pal[4] = [0.2154516905067, 0.295, 0.5078367867510504]
pal[5] = [0.78, 0.51, 0.1892045842]
pal[6] = pal[7] = pal[8] = pal[9] = [0.12071162840208301, 0.4357915995132193, 0.2463679091477368]
pal[10] = pal[11] = pal[12] = pal[13] = [0.7, 0.7, 0.65]

tiles = sns.FacetGrid(df, row="tiles", hue="tiles", aspect=15, height=.5, palette=pal)

# densities
tiles.map(sns.kdeplot, "timescale",
      bw_adjust=.5, clip_on=False,
      fill=True, alpha=1, linewidth=1.5)
tiles.map(sns.kdeplot, "timescale", clip_on=False, color="w", lw=2, bw_adjust=.5)
tiles.refline(y=0, linewidth=2, linestyle="-", color=None, clip_on=False)


def label(x, color, label):
    ax = plt.gca()
    ax.text(0, .2, label, fontweight="bold", color=color,
            ha="left", va="center", transform=ax.transAxes)
    
tiles.map(label, "timescale")

# overlapping tiles
tiles.figure.subplots_adjust(hspace=-0.15)

# formatting for visibility 
tiles.set_titles("")
tiles.despine(bottom=True, left=True)
plt.title('4-Cycles without OH', y=13.5, fontsize=14)

ticks = plt.xticks()
xtick_string = [rf'$10^{str({i})}$' for i in range(-10, 25, 5)]

tiles.set(xticks=range(-10, 25, 5))
tiles.set_xticklabels(xtick_string)


tiles.set(xlabel='Timescale (seconds/cycle)')
tiles.set(yticks=[], ylabel="")

plt.show
