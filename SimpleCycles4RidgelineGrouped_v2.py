#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 28 02:13:45 2023

@author: emywli
"""

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})

# Pulling the data
# cycles4data = open('/Users/emywli/Desktop/Research/TwilightAmazon_L10_4cyclestimes.csv','r')
# data4 = cycles4data.read()
# lines4 = data4.split('\n')
# cycles4data.close()

# lines4.pop(0)
# # lines4.pop(-1)

# cycleseight = []
# for i in range(0, len(lines8)): 
#     temp = lines8[i].split(',')
#     cycleseight.append(temp)

# getting the timescales from stored sorted cycles .csv files
cycles1data = open('/Users/emywli/Desktop/Research/L1 Noon Cycles/Amazon_L1_4cyclestimes_simple.csv','r')
data1 = cycles1data.read()
lines1 = data1.split('\n')
cycles1data.close()

lines1.pop(0)
lines1.pop(-1)

times1 = []
for i in range(0, len(lines1)): 
    temp = lines1[i].split(',"')
    # print(temp)
    times1.append(float(temp[0]))
    
cycles2data = open('/Users/emywli/Desktop/Research/L1 Noon Cycles/AtlanticOcean_L1_4cyclestimes_simple.csv','r')
data2 = cycles2data.read()
lines2 = data2.split('\n')
cycles2data.close()

lines2.pop(0)
lines2.pop(-1)

times2 = []
for i in range(0, len(lines2)): 
    temp = lines2[i].split(',"')
    # print(temp)
    times2.append(float(temp[0]))

cycles3data = open('/Users/emywli/Desktop/Research/L1 Noon Cycles/Beijing_L1_4cyclestimes_simple.csv','r')
data3 = cycles3data.read()
lines3 = data3.split('\n')
cycles3data.close()

lines3.pop(0)
lines3.pop(-1)

times3 = []
for i in range(0, len(lines3)): 
    temp = lines3[i].split(',"')
    # print(temp)
    times3.append(float(temp[0]))
    
cycles4data = open('/Users/emywli/Desktop/Research/L1 Noon Cycles/Borneo_L1_4cyclestimes_simple.csv','r')
data4 = cycles4data.read()
lines4 = data4.split('\n')
cycles4data.close()

lines4.pop(0)
lines4.pop(-1)

times4 = []
for i in range(0, len(lines4)): 
    temp = lines4[i].split(',"')
    # print(temp)
    times4.append(float(temp[0]))

cycles5data = open('/Users/emywli/Desktop/Research/L1 Noon Cycles/CapeGrim_L1_4cyclestimes_simple.csv','r')
data5 = cycles5data.read()
lines5 = data5.split('\n')
cycles5data.close()

lines5.pop(0)
lines5.pop(-1)

times5 = []
for i in range(0, len(lines5)): 
    temp = lines5[i].split(',"')
    # print(temp)
    times5.append(float(temp[0]))
    
cycles6data = open('/Users/emywli/Desktop/Research/L1 Noon Cycles/Congo_L1_4cyclestimes_simple.csv','r')
data6 = cycles6data.read()
lines6 = data6.split('\n')
cycles6data.close()

lines6.pop(0)
lines6.pop(-1)

times6 = []
for i in range(0, len(lines6)): 
    temp = lines6[i].split(',"')
    # print(temp)
    times6.append(float(temp[0]))
    
cycles7data = open('/Users/emywli/Desktop/Research/L1 Noon Cycles/ElDjouf_L1_4cyclestimes_simple.csv','r')
data7 = cycles7data.read()
lines7 = data7.split('\n')
cycles7data.close()

lines7.pop(0)
lines7.pop(-1)

times7 = []
for i in range(0, len(lines7)): 
    temp = lines7[i].split(',"')
    # print(temp)
    times7.append(float(temp[0]))
    
cycles8data = open('/Users/emywli/Desktop/Research/L1 Noon Cycles/Graciosa_L1_4cyclestimes_simple.csv','r')
data8 = cycles8data.read()
lines8 = data8.split('\n')
cycles8data.close()

lines8.pop(0)
lines8.pop(-1)

times8 = []
for i in range(0, len(lines8)): 
    temp = lines8[i].split(',"')
    # print(temp)
    times8.append(float(temp[0]))
    
cycles9data = open('/Users/emywli/Desktop/Research/L1 Noon Cycles/IndianOcean_L1_4cyclestimes_simple.csv','r')
data9 = cycles9data.read()
lines9 = data9.split('\n')
cycles9data.close()

lines9.pop(0)
lines9.pop(-1)

times9 = []
for i in range(0, len(lines9)): 
    temp = lines9[i].split(',"')
    # print(temp)
    times9.append(float(temp[0]))
    
cycles10data = open('/Users/emywli/Desktop/Research/L1 Noon Cycles/Kinshasa_L1_4cyclestimes_simple.csv','r')
data10 = cycles10data.read()
lines10 = data10.split('\n')
cycles10data.close()

lines10.pop(0)
lines10.pop(-1)

times10 = []
for i in range(0, len(lines10)): 
    temp = lines10[i].split(',"')
    # print(temp)
    times10.append(float(temp[0]))
    
cycles11data = open('/Users/emywli/Desktop/Research/L1 Noon Cycles/LosAngeles_L1_4cyclestimes_simple.csv','r')
data11 = cycles11data.read()
lines11 = data11.split('\n')
cycles11data.close()

lines11.pop(0)
lines11.pop(-1)

times11 = []
for i in range(0, len(lines11)): 
    temp = lines11[i].split(',"')
    # print(temp)
    times11.append(float(temp[0]))
    
cycles12data = open('/Users/emywli/Desktop/Research/L1 Noon Cycles/McMurdo_L1_4cyclestimes_simple.csv','r')
data12 = cycles12data.read()
lines12 = data12.split('\n')
cycles12data.close()

lines12.pop(0)
lines12.pop(-1)

times12 = []
for i in range(0, len(lines12)): 
    temp = lines12[i].split(',"')
    # print(temp)
    times12.append(float(temp[0]))
    
cycles13data = open('/Users/emywli/Desktop/Research/L1 Noon Cycles/Ozarks_L1_4cyclestimes_simple.csv','r')
data13 = cycles13data.read()
lines13 = data13.split('\n')
cycles13data.close()

lines13.pop(0)
lines13.pop(-1)

times13 = []
for i in range(0, len(lines13)): 
    temp = lines13[i].split(',"')
    # print(temp)
    times13.append(float(temp[0]))
    
cycles14data = open('/Users/emywli/Desktop/Research/L1 Noon Cycles/PacificOcean_L1_4cyclestimes_simple.csv','r')
data14 = cycles14data.read()
lines14 = data14.split('\n')
cycles14data.close()

lines14.pop(0)
lines14.pop(-1)

times14 = []
for i in range(0, len(lines14)): 
    temp = lines14[i].split(',"')
    # print(temp)
    times14.append(float(temp[0]))
    
cycles15data = open('/Users/emywli/Desktop/Research/L1 Noon Cycles/Paris_L1_4cyclestimes_simple.csv','r')
data15 = cycles15data.read()
lines15 = data15.split('\n')
cycles15data.close()

lines15.pop(0)
lines15.pop(-1)

times15 = []
for i in range(0, len(lines15)): 
    temp = lines15[i].split(',"')
    # print(temp)
    times15.append(float(temp[0]))
    
cycles16data = open('/Users/emywli/Desktop/Research/L1 Noon Cycles/Utqiagvik_L1_4cyclestimes_simple.csv','r')
data16 = cycles16data.read()
lines16 = data16.split('\n')
cycles16data.close()

lines16.pop(0)
lines16.pop(-1)

times16 = []
for i in range(0, len(lines16)): 
    temp = lines16[i].split(',"')
    # print(temp)
    times16.append(float(temp[0]))
    

# times4arr = np.array(times4)
# times6arr = np.array(times6)
# times8arr = np.array(times8)
# x = np.array(times4 + times6 + times8)
# tiles = []
# temp = 4
# while len(tiles) < 3:
#     tiles.append(str(temp))
#     temp += 2
    
tiles = []
# tiles = ['4', '6', '8']
for i in range(0, len(times2)): 
    tiles.append('Atlantic Ocean')   
for i in range(0, len(times9)): 
    tiles.append('Indian Ocean')  
for i in range(0, len(times14)): 
    tiles.append('Pacific Ocean')  
for i in range(0, len(times8)): 
    tiles.append('Graciosa')  
for i in range(0, len(times5)): 
    tiles.append('Kennaook')  
for i in range(0, len(times7)): 
    tiles.append('El Djouf')  
for i in range(0, len(times1)): 
    tiles.append('Amazon')
for i in range(0, len(times4)): 
    tiles.append('Borneo')  
for i in range(0, len(times6)): 
    tiles.append('Congo')  
for i in range(0, len(times13)): 
    tiles.append('Ozarks')  
for i in range(0, len(times3)): 
    tiles.append('Beijing')  
for i in range(0, len(times10)): 
    tiles.append('Kinshasa')  
for i in range(0, len(times11)): 
    tiles.append('Los Angeles')   
for i in range(0, len(times15)): 
    tiles.append('Paris')  
for i in range(0, len(times16)): 
    tiles.append('Utqiagvik')     
for i in range(0, len(times12)): 
    tiles.append('McMurdo Station') 

timescalelin = times2+times9+times14+times8+times5+times7+times1+times4+times6+times13+times3+times10+times11+times15+times16+times12
timescalelin = np.array(timescalelin)
# timescale = timescalelin
timescale = np.log10(timescalelin)
# g = np.tile(list("ABCDEFGHIJ"), 50)
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

# plt.savefig('/Users/emywli/Desktop/Research/RidgelineLocationsL1NoonGroupedv2.pdf', bbox_inches="tight")
