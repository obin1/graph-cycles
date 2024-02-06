#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 28 02:13:45 2023

@author: emywli
"""
# https://matplotlib.org/stable/gallery/text_labels_and_annotations/date.html#sphx-glr-gallery-text-labels-and-annotations-date-py

import numpy as np
import networkx as nx
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from pathlib import Path
from datetime import datetime, timedelta
import decimal
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

spc_list = species

# targetcycle = {'NO2', 'R760', 'NO', 'R13'}
# for i in range(0, len(cyclesfour)):
#     if set(cyclesfour[i]) == targetcycle:
#         targetcycle = list(cyclesfour[i])


# list to store the times of day 
localtime = []

# initializing dictionary to hold all the different locations and the time series for each 
plotcyctimes = {
   'Primary Photolytic': ([],[],[]),
   'HO2 Oxidation': ([],[],[]),
   'IHOO1': ([],[],[]),
   'IHOO4': ([],[],[])
}

# list of cycles, based on participating reactions, to compare 
cycleslist = [[760, 200, 13], [22, 760], [420, 760], [422, 760], [424, 760], [426, 760]]

directory = 'TimeSeries_Samples'
loc = Path(directory)
locfolders = [x for x in loc.iterdir() if x.is_dir()]

for folder in locfolders: 
    location = (folder.name).split('_')
    locname = location[0]
    if locname == 'Kinshasa':
        locfiles = Path(folder).rglob('*.txt')
        # determining how many hours to add/subtract from the file GMT time to obtain local time, Graciosa and El Djouf already on GMT
        hoursdelta = timedelta(hours=0)
        if locname == 'Atlantic Ocean': 
            hoursdelta = timedelta(hours = -3)
        elif locname == 'Indian Ocean': 
            hoursdelta = timedelta(hours = 6)
        elif locname == 'Pacific Ocean': 
            hoursdelta = timedelta(hours = -7)
        elif locname == 'Kennaook': 
            hoursdelta = timedelta(hours = 10)
        elif locname == 'Amazon': 
            hoursdelta = timedelta(hours = -4)
        elif locname == 'Borneo': 
            hoursdelta = timedelta(hours = 8)
        elif locname == 'Congo' :
            hoursdelta = timedelta(hours = 1)
        elif locname == 'Ozarks': 
            hoursdelta = timedelta(hours = -5)
        elif locname == 'Beijing': 
            hoursdelta = timedelta(hours = 8)
        elif locname == 'Kinshasa': 
            hoursdelta = timedelta(hours = 1)
        elif locname == 'Los Angeles': 
            hoursdelta = timedelta(hours = -7)
        elif locname == 'Paris':
            hoursdelta = timedelta(hours = 2)
        elif locname == 'Utqiagvik':
            hoursdelta = timedelta(hours = -8)
        elif locname == 'McMurdo Station': 
            hoursdelta == timedelta(hours = 12)
        # lists to store the timescale data of each target cycles for every 15-minute interval of the day
        cycle1times = []
        cycle2times = []
        cycle3times = []
        cycle4times = []
        timescalelin = []
        for file in locfiles:
            full_name = Path(file).name
            full_name = full_name.split('_')
            newname = full_name[3].split('.')
            if 'Twilight' in locname: 
                locname = locname.replace('Twilight', '')
            timestampstr = newname[0]
            # converting GMT time to local time 
            time_object = datetime.strptime(timestampstr, '%H%M')
            time_object2 = time_object + hoursdelta
            timenew = time_object2.time()
            loctime = timenew.strftime('%H:%M')
            sample_file = open(file, 'rt')
            sampledata = sample_file.read()
            sample = sampledata.split('\n ')
            sample_file.close()
            # list to hold the timescales for all cycles in cycleslist for each file (time of day) 
            timescalesraw = [] 
            for i in range(0, len(cycleslist)): 
                cyclet = 0
                for r in range(0, len(sample)): 
                    if sample[r][0] == 'A': 
                        temp = sample[r].split()
                        if int(temp[1]) in cycleslist[i]: 
                            if float(temp[4]) == float(0):
                                cyclet += float('inf')
                            else:
                                cyclet += 1/(float(temp[4]))
                                # print(cyclet)
                timescalesraw.append(cyclet) 
            for t in range(0, len(timescalesraw)):
                # variables for determining characteristic timescales of sets of parallel/simultaneous cycles 
                cycrate_a = 0
                cycrate_b = 0
                combinedr = 0 
                parallelt = 0
                if t == 0 or t == 1: 
                    varname = 'cycle' + str(t+1) +'times'
                    globals()[varname].append([loctime, timescalesraw[t]]) 
                elif t == 2: 
                    cycrate_a = 1/(float(timescalesraw[t]))
                    cycrate_b = 1/(float(timescalesraw[t+1])) 
                    combinedr = cycrate_a + cycrate_b
                    if combinedr != float(0):
                        parallelt = 1/(float(combinedr)) 
                    else: 
                        parallelt = float('inf')
                    cycle3times.append([loctime, parallelt]) 
                elif t == 4: 
                    cycrate_a = 1/(float(timescalesraw[t]))
                    cycrate_b = 1/(float(timescalesraw[t+1])) 
                    combinedr = cycrate_a + cycrate_b
                    if combinedr != float(0):
                        parallelt = 1/(float(combinedr)) 
                    else: 
                        parallelt = float('inf')
                    cycle4times.append([loctime, parallelt]) 
            cycle1sorted = sorted(cycle1times, key=lambda x: x[0])
            cycle2sorted = sorted(cycle2times, key=lambda x: x[0])
            cycle3sorted = sorted(cycle3times, key=lambda x: x[0])
            cycle4sorted = sorted(cycle4times, key=lambda x: x[0])
        for s in range(1, 5): 
            cyclevals = 'cycle' + str(s) +'sorted'
            list1 = []
            list2 = []
            for j in range(0, len(globals()[cyclevals])):
                list1.append(globals()[cyclevals][j][0])
                list2.append(globals()[cyclevals][j][1])
            list2arr = np.array(list2)
            list2new = np.log10(list2arr)
            if s == 1: 
                cycname = 'Primary Photolytic'
                plotcyctimes[cycname] = (list1, list2new, [0.78, 0.51, 0.1892045842])
            elif s == 2: 
                cycname = 'HO2 Oxidation'
                plotcyctimes[cycname] = (list1, list2new, [0.2154516905067, 0.295, 0.5078367867510504])
            elif s == 3: 
                cycname = 'IHOO1'
                plotcyctimes[cycname] = (list1, list2new, [0.12071162840208301, 0.4357915995132193, 0.2463679091477368])
            elif s == 4: 
                cycname = 'IHOO4'
                plotcyctimes[cycname] = (list1, list2new, [0.10071162840208301, 0.63995132193, 0.20679091477368])

for label, (x, y, c) in plotcyctimes.items():
    plt.plot(x, y, label = label, color=c, linewidth=0.8)
# Adding legend, x and y labels, and title for the lines
plt.legend(facecolor='white', framealpha=0, fontsize=8, bbox_to_anchor=(1, 1))
plt.xlabel("Local Time")
plt.xticks(fontsize=4, rotation=90)
# plt.tick_params(which='major', direction='out', length=5, width=10, color='black')
plt.ylabel('Cycle Timescale (log(seconds/(molec/cc)))')
plt.title('Comparison of Cycle Timescales in Kinshasa over a 24-Hour Period')
# Displaying the plot 
# plt.show
plt.savefig('/Users/emywli/Desktop/Research/Time Series Cycle Comparisons/TimeSeriesL1July_CycleComparison_Kinshasa.png', dpi=300, bbox_inches="tight")



