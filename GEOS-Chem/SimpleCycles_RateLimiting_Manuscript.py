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
import matplotlib.colors as colors
from matplotlib.pyplot import figure
sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})


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

# Create DiGraph
B = nx.DiGraph()

cyc4_horiz = np.array([])
cyc4_vert = np.array([])

cyc6_horiz = np.array([])
cyc6_vert = np.array([])

cyc8_horiz = np.array([])
cyc8_vert = np.array([])

# starting big loop, each iteration to get the plot inputs for different cycle lengths 
for l in range(4, 10, 2): 
    print(l)
    cyclength = l  
    thisfile = 'GCv14_' + str(cyclength) + 'cycles_simple.csv'
    # Reading in stored cycles from csv 
    cycles = open(thisfile,'r')
    datacyc = cycles.read()
    linescyc = datacyc.split('\n')
    cycles.close()
    # removing first and last lines which are blank in the saved CSVs
    linescyc.pop(0)
    linescyc.pop(-1)
    cycleslist = []
    for i in range(0, len(linescyc)): 
        temp = linescyc[i].split(',')
        cycleslist.append(temp)
        
    totalcycles = cycleslist
    
    # list to store the data for cycle timescale distributions
    timescalelin = []
    
    # list to store max timescales in each cycle, across all 16 locations as well
    maxtimescalelin = []
    
    # list to store the repeat location names corresponding to the number of timescale values for that location
    tiles = []
    
    # list to store all the longest-to-full cycle timescale ratios 
    maxratios = []
    
    # location names as they are in the file names (all at noon local time for this analysis case), in order of how we want them to show up on graph
    locations = [['AtlanticOcean', '1500'], ['IndianOcean', '0600'], ['PacificOcean', '1900'], ['Graciosa', '1200'], ['CapeGrim', '0200'], ['ElDjouf', '1200'], ['Amazon', '1600'], ['Borneo', '0400'], ['Congo', '1100'], ['Ozarks', '1700'], ['Beijing', '0400'], ['Kinshasa', '1100'], ['LosAngeles', '1900'], ['Paris', '1000'], ['Utqiagvik', '2000'], ['McMurdo', '0000']]
    # location names as how we want them to appear on the graph
    graphlocations = ['Atlantic Ocean', 'Indian Ocean', 'Pacific Ocean', 'Graciosa', 'Kennaook', 'El Djouf', 'Amazon', 'Borneo', 'Congo', 'Ozarks', 'Beijing', 'Kinshasa', 'Los Angeles', 'Paris', 'Utqiagvik', 'McMurdo Station']
    # filename = '/Users/emywli/Desktop/Research/L1 Noon Samples/AtlanticOcean_L1_20180702_2000.txt'
    for location in locations: 
        filename = '/Users/emywli/Desktop/Research/graph-cycles/GEOS-Chem/L1 Noon Samples/' + str(location[0]) + '_L1_20180702_' + str(location[1]) + '.txt'
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
        # list to hold residence times for each cycle in totalcycles
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
            # add time per reaction to get residence time for the full cycle on current iteration 
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
            if cyclest[i][0] != float('inf'): 
                cyclest_noninf.append([cyclest[i][0], cyclest[i][1], cyclest[i][2], cyclest[i][3]])
        cyclest_noninf.sort()
        ratelimitingtimescales = []
        fulltimescales = []
        ratios = []
        for i in range(0, len(cyclest_noninf)): 
            ratelimitingtimescales.append(float(cyclest_noninf[i][1]))
            fulltimescales.append(float(cyclest_noninf[i][0]))
            # keeping track of the cycle to which each rate limiting ratio belongs, and in which location
            ratios.append([float(cyclest_noninf[i][2]), cyclest_noninf[i][3], locations.index(location)])
        timescalelin += fulltimescales
        maxtimescalelin += ratelimitingtimescales
        maxratios += ratios
    
    maxratiosvals = []
    for i in range(0, len(maxratios)):
        maxratiosvals.append(maxratios[i][0])
        
    maxtimescalelin = np.array(maxtimescalelin)
    timescalelin = np.array(timescalelin)
    maxratiosarr = np.array(maxratiosvals)
    hname = 'cyc' + str(l) + '_horiz'
    vname = 'cyc' + str(l) + '_vert'
    globals()[hname] = np.log10(timescalelin)
    globals()[vname] = maxratiosarr
    print(len(maxratiosarr))

# generating the figure of 3 subfigures 
fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, figsize=(18, 6), sharey=True, layout='constrained')

ax1.hexbin(cyc4_horiz, cyc4_vert, gridsize = [57,16], cmap='viridis_r', mincnt=1, linewidths=0.1, norm=colors.LogNorm()) 
ax1.set_xlabel("Full Cycle Timescale (seconds/cycle)", fontsize=21, labelpad=9)
ax1.set_title('Two-Reaction Cycles', fontsize=26) 
ax1.set_ylim(0.0, 1.05)
xtick_string = [rf'$10^{str({i})}$' for i in range(-10, 45, 10)]
ax1.set_xticks(range(-10, 45, 10), labels = xtick_string, fontsize=21)
ytx = ax1.get_yticks()
ytxstr = [f"{t:.1f}" for t in ytx]
ax1.set_yticks(ytx, labels=ytxstr, fontsize=21)
ax1.tick_params(direction='out', color='black', length=5, bottom=True, left=True)
ax1.set_ylabel("Rate-Limiting Reaction Timescale/\nFull Timescale Ratio", fontsize=20)
ax1.axvline(x = 0, color = 'gray', label = 'axvline - full height', linewidth=0.5)
ax1.set(xlim=(-12, 45))

ax2.hexbin(cyc6_horiz, cyc6_vert, gridsize = [48,16], cmap='plasma_r', mincnt=1, linewidths=0.1, norm=colors.LogNorm())
ax2.set_xlabel("Full Cycle Timescale (seconds/cycle)", fontsize=21, labelpad=9)
ax2.set_title('Three-Reaction Cycles', fontsize=26) 
ax2.set_ylim(0.0, 1.05)
xtick_string = [rf'$10^{str({i})}$' for i in range(-10, 45, 10)]
ax2.set_xticks(range(-10, 45, 10), labels = xtick_string, fontsize=21)
ax2.tick_params(direction='out', color='black', length=5, bottom=True)
ax2.set_xticklabels(xtick_string)
ax2.axvline(x = 0, color = 'gray', label = 'axvline - full height', linewidth=0.5)
ax2.set(xlim=(-12, 45))

ax3.hexbin(cyc8_horiz, cyc8_vert, gridsize = [45,15], cmap='PuBuGn', mincnt=1, linewidths=0.1, norm=colors.LogNorm())
ax3.set_xlabel("Full Cycle Timescale (seconds/cycle)", fontsize=21, labelpad=9)
ax3.set_title('Four-Reaction Cycles', fontsize=23) 
ax3.set_ylim(0.0, 1.05)
xtick_string = [rf'$10^{str({i})}$' for i in range(-10, 45, 10)]
ax3.set_xticks(range(-10, 45, 10), labels = xtick_string, fontsize=21)
ax3.tick_params(direction='out', color='black', length=5, bottom=True)
ax3.set_xticklabels(xtick_string)
ax3.axvline(x = 0, color = 'gray', label = 'axvline - full height', linewidth=0.5)
ax3.set(xlim=(-12, 45))

# setting the colorbars for the subfigures 
cbar = fig.colorbar(ax1.collections[0], ax=[ax1], orientation='horizontal', pad=0.08)
cbar.ax.minorticks_off()
cbar.ax.tick_params(labelsize=21)

cbar = fig.colorbar(ax2.collections[0], ax=[ax2], orientation='horizontal', pad=0.08)
cbar.ax.tick_params(labelsize=21)
cbar.ax.minorticks_off()  

cbar = fig.colorbar(ax3.collections[0], ax=[ax3], orientation='horizontal', pad=0.08)
cbar.ax.tick_params(labelsize=21)
cbar.ax.minorticks_off()

plt.show
# plt.savefig('RateLimitedCycles_HexBinFigure.png', dpi=300, bbox_inches="tight")
