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
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import matplotlib.patches as mpatches
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

reactions = ['R' + str(i) for i in range(1,914)]

# Create DiGraph
B = nx.DiGraph()


# Reading in stored cycles from csv 
cycles = open('GCv14_4cycles_simple.csv','r')
data = cycles.read()
lines = data.split('\n')
cycles.close()

lines.pop(0)
lines.pop(-1)

cyclesfour = []
for i in range(0, len(lines)): 
    temp = lines[i].split(',')
    cyclesfour.append(temp)
    
totalcycles = cyclesfour

# list storing lists of future timescalelin elements, for location ordering purposes 
timescalelinlist = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]

# list to store the data for distributions
timescalelin = []

# list to store lists of future tiles elements, for location ordering purposes
tilesog = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]

# list to store the repeat location names corresponding to the number of timescale values for that location
tiles = []

# list to hold the sum of all inverse timescales (parallel rate of cycles) and corresponding location 
# [location, time, sum of all cycle rates]
sumsloc = []

loc0vals = []
loc1vals = []
loc2vals = []
loc3vals = []
loc4vals = []
loc5vals = []
loc6vals = []
loc7vals = []
loc8vals = []
loc9vals = []
loc10vals = []
loc11vals = []
loc12vals = []
loc13vals = []
loc14vals = []
loc15vals = []
loc16vals = []

# location names as they are in the file names, in order of how we want them to show up on graph
# locations = [['AtlanticOcean', '1500'], ['IndianOcean', '0600'], ['PacificOcean', '1900'], ['Graciosa', '1200'], ['CapeGrim', '0200'], ['ElDjouf', '1200'], ['Amazon', '1600'], ['Borneo', '0400'], ['Congo', '1100'], ['Ozarks', '1700'], ['Beijing', '0400'], ['Kinshasa', '1100'], ['LosAngeles', '1900'], ['Paris', '1000'], ['Utqiagvik', '2000'], ['McMurdo', '0000']]
locations = ['AtlanticOcean', 'IndianOcean', 'PacificOcean', 'Graciosa', 'CapeGrim', 'ElDjouf', 'Amazon', 'Borneo', 'Congo', 'Ozarks', 'Beijing', 'Kinshasa', 'LosAngeles', 'Paris', 'Utqiagvik', 'McMurdo']
# location names as how we want them to appear on the graph
graphlocations = ['Atlantic Ocean', 'Indian Ocean', 'Pacific Ocean', 'Graciosa', 'Kennaook', 'El Djouf', 'Amazon', 'Borneo', 'Congo', 'Ozarks', 'Beijing', 'Kinshasa', 'Los Angeles', 'Paris', 'Utqiagvik', 'McMurdo Station']
# 16 colors based on location classification
colorslist = [[0.2154516905067, 0.295, 0.5078367867510504], [0.2154516905067, 0.295, 0.5078367867510504], [0.2154516905067, 0.295, 0.5078367867510504], [0.2154516905067, 0.295, 0.5078367867510504], [0.2154516905067, 0.295, 0.5078367867510504], [0.78, 0.51, 0.1892045842], [0.12071162840208301, 0.4357915995132193, 0.2463679091477368], [0.12071162840208301, 0.4357915995132193, 0.2463679091477368], [0.12071162840208301, 0.4357915995132193, 0.2463679091477368], [0.12071162840208301, 0.4357915995132193, 0.2463679091477368], [0.7, 0.7, 0.65], [0.7, 0.7, 0.65], [0.7, 0.7, 0.65], [0.7, 0.7, 0.65], [0.5184891135337089, 0.7194769933793438, 0.7518803726887796], [0.5632111255041908, 0.758620966612444, 0.7764634182455044]]
# 16 markers corresponding to the 16 locations for the repeat location classifications 
markerslist = ['o', '^', 's', 'p', 'd', 'o', 'o', '^', 's', 'p', 'o', '^', 's', 'p', 'o', '^']

directory = 'Samples by Hour'
loc = Path(directory)
locfolders = [x for x in loc.iterdir() if x.is_dir()]

# each folder holding 16 files for the 16 different locations at that time 
for folder in locfolders: 
    folname = folder.name.split('_')
    timestampstr = folname[1]
    time_object = datetime.strptime(timestampstr, '%H%M')
    timenew = time_object.time()
    loctime = timenew.strftime('%H:%M')
    print(folder.name)
    locfiles = Path(folder).rglob('*.txt')
    for file in locfiles: 
        full_name = Path(file).name
        location = full_name.split('_')
        locname = location[0]
        if 'Twilight' in locname: 
            locname = locname.replace('Twilight', '')
        # print(full_name)
        # print(locname)
        sample_file = open(file, "rt")
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
        # list to hold timescale for each cycle in totalcycles
        cyclest = []
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
            # add time per reaction to get timescale for the full cycle on current iteration 
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
        # removing the timescale that are 'inf' for reasonability, noting that inf would skew it higher
        cyclest_noninf = []
        for i in range(0, len(cyclest)):
            if np.log10(cyclest[i][0]) <= 20:
                cyclest_noninf.append([cyclest[i][0], cyclest[i][1], cyclest[i][2], cyclest[i][3]])
        cyclest_noninf.sort()
        times = []
        for i in range(0, len(cyclest_noninf)): 
            times.append(float(cyclest_noninf[i][0]))
        sumrates = 0.0
        # sumratesperc variable for summing only the inverse timescales that fall within a certain range (e.g. 25th to 75th percentile of timescales)
        sumratesperc = 0.0
        sumlist = []
        for j in range(0, len(times)): 
            sumlist.append(float(1/times[j]))
            sumrates += float(1/times[j])
        # below percentile numbers can be modified to examine different portions of the timescale distribution -- and how much the time series changes depending on what percentile we look at closely across locations... 
        # ... i.e. focusing on the slowest or fastest timescales in the different locations changing the differences in their total number of cycles
        perc25 = np.percentile(sumlist, 25)
        perc25_idx = np.argmin(sumlist>perc25)
        perc75 = np.percentile(sumlist, 75)
        perc75_idx = np.argmax(sumlist<perc75)
        sumlistnew = []
        for t in range(perc75_idx, perc25_idx): 
            sumlistnew.append(sumlist[t])
            sumratesperc += sumlist[t]
        # getting index of current iteration's location to put into tiles list
        index = locations.index(locname)
        locofficial = graphlocations[index]
        varname = 'loc' + str(index) +'vals'
        globals()[varname].append([loctime, sumratesperc, sumrates, perc75, perc75_idx, perc25, perc25_idx]) 
        # print(index)
        # tiles += len(times) * [graphlocations[index]]
        # timescalelinlist[index] = times
        # tilesog[index] = (len(times) * [graphlocations[index]])

plotloctimes = {
    'Atlantic Ocean': ([],[],[],[]),
    'Indian Ocean': ([],[],[],[]),
    'Pacific Ocean': ([],[],[],[]),
    'Graciosa': ([],[],[],[]),
    'Kennaook': ([],[],[],[]),
    'El Djouf': ([],[],[],[]),
    'Amazon': ([],[],[],[]),
    'Borneo': ([],[],[],[]),
    'Congo': ([],[],[],[]),
    'Ozarks': ([],[],[],[]),
    'Beijing': ([],[],[],[]),
    'Kinshasa': ([],[],[],[]),
    'Los Angeles': ([],[],[],[]),
    'Paris': ([],[],[],[]),
    'Utqiagvik': ([],[],[],[]),
    'McMurdo Station': ([],[],[],[]),
}

loc0sorted = sorted(loc0vals, key=lambda x: x[0])
loc1sorted = sorted(loc1vals, key=lambda x: x[0])
loc2sorted = sorted(loc2vals, key=lambda x: x[0])
loc3sorted = sorted(loc3vals, key=lambda x: x[0])
loc4sorted = sorted(loc4vals, key=lambda x: x[0])
loc5sorted = sorted(loc5vals, key=lambda x: x[0])
loc6sorted = sorted(loc6vals, key=lambda x: x[0])
loc7sorted = sorted(loc7vals, key=lambda x: x[0])
loc8sorted = sorted(loc8vals, key=lambda x: x[0])
loc9sorted = sorted(loc9vals, key=lambda x: x[0])
loc10sorted = sorted(loc10vals, key=lambda x: x[0])
loc11sorted = sorted(loc11vals, key=lambda x: x[0])
loc12sorted = sorted(loc12vals, key=lambda x: x[0])
loc13sorted = sorted(loc13vals, key=lambda x: x[0])
loc14sorted = sorted(loc14vals, key=lambda x: x[0])
loc15sorted = sorted(loc15vals, key=lambda x: x[0])

timesofday = []
for i in range(0, len(graphlocations)):
    locgraph = graphlocations[i]
    sortedname = 'loc' + str(i) + 'sorted'
    timesofday = []
    # for the total number of cycles occurring at each hour 
    numcyclesoccurring = []
    for j in range(0, 24): 
        timesofday.append(globals()[sortedname][j][0])
        numcyclesoccurring.append(globals()[sortedname][j][2])
    plotloctimes[locgraph] = (timesofday, numcyclesoccurring, colorslist[i], markerslist[i])


plt.figure(figsize=(15,8))
# for shading in the regions for each category bounded by the minimum sumrate for that category and the maximum, at each given time
# locations (by corresponding indices) in their categories
groups = [[0,1,2,3,4],[5],[6,7,8,9],[10,11,12,13],[14],[15]]
fillcolors = [[0.2154516905067, 0.295, 0.5078367867510504],[0.78, 0.51, 0.1892045842], [0.12071162840208301, 0.4357915995132193, 0.2463679091477368],[0.7, 0.7, 0.65],[0.5184891135337089, 0.7194769933793438, 0.7518803726887796],[0.7332111255041908, 0.58620966612444, 0.7864634182455044]]
loclabels = ['Ocean', 'Desert', 'Biogenic', 'Urban', 'Arctic', 'Polar']
# patches = []
for i in range(0, len(groups)):
    ymins = []
    ymaxs = []
    groupavgs = []
    for j in range(0, 24): 
        tempcompare = []
        sumcompare = 0
        avg = 0 
        for k in range(0, len(groups[i])):
            sortedname = 'loc' + str(groups[i][k]) + 'sorted'
            tempcompare.append(globals()[sortedname][j][2]) 
        # print(tempcompare)
        for m in range(0, len(tempcompare)): 
            sumcompare += tempcompare[m]
        avg = sumcompare / len(tempcompare)
        groupavgs.append(avg)
        tempmin = min(tempcompare)
        tempmax = max(tempcompare)
        ymins.append(tempmin)
        ymaxs.append(tempmax)
    # print(ymins)
    # print(ymaxs)
    if i == 0 or i == 2 or i == 3:
        plt.plot(timesofday, groupavgs, label = loclabels[i], linewidth=3 , color=fillcolors[i])
        plt.fill_between(timesofday, ymins, ymaxs, color=fillcolors[i], alpha=0.25, linewidth=0)
plt.legend(facecolor='white', framealpha=0, fontsize=20, bbox_to_anchor=(1.2, 1))
plt.yscale("log")
plt.xlabel("Local Time", fontsize=22, labelpad=15)
plt.xticks(fontsize=18, rotation=90)
plt.yticks(fontsize=22)
plt.ylabel('Total Cycles Rate (cycles/sec)', fontsize=22)        
        
        
plt.show
# plt.savefig('TimeSeriesL1_Figure.png', dpi=300, bbox_inches="tight")


# testing, for each location, how much a certain percentile range covers out of all cycles occurring at each time
# coverageratios = []
# for i in range(0, len(graphlocations)):
#     tempratios = []
#     locgraph = graphlocations[i]
#     sortedname = 'loc' + str(i) + 'sorted'
#     for j in range(0, len(globals()[sortedname])): 
#         temp = float(globals()[sortedname][j][1])/float(globals()[sortedname][j][2])
#         tempratios.append(temp)
#     coverageratios.append(tempratios)
      
# avgratios = []
# for i in range(0, len(coverageratios)): 
#     totalratios = 0.0
#     for j in range(0, len(coverageratios[i])): 
#         totalratios+= coverageratios[i][j]
#     avgratio = totalratios/24
#     avgratios.append(avgratio)