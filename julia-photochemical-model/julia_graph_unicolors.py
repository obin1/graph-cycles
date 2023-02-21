#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 11 12:39:10 2023

@author: psturm
"""
import numpy as np
import networkx as nx
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
import random as rd
# Note bipartite tools are not imported into the networkx namespace
from networkx.algorithms import bipartite


# Reference model 
# from Sturm and Wexler (2022): https://doi.org/10.5194/gmd-15-3417-2022
#  Species
#         1) O3
#         2) NO
#         3) NO2
#         4) HCHO
#         5) HO2
#         6) HO2H
#         steady state species: .HO, O
#         buildup species: HNO3, CO, H2

#  Reactions                         
#         1) NO2 + HV = NO + O                             
#         2) O + O2 = O3
#         3) O3 + NO = NO2 + O2
#         4) HCHO + HV = 2 HO2. + CO
#         5) HCHO + HV = H2 + CO
#         6) HCHO + HO. = HO2. + CO + H2O
#         7) HO2. + NO = HO. + NO2
#         8) HO. + NO2 = HNO3
#         9) HO2H + HV = 2 HO.
#        10) HO2H + HO. = H2O + HO2.



A = np.matrix('0, 1,-1, 0, 0, 0, 0, 0, 0, 0;\
               1, 0,-1, 0, 0, 0,-1, 0, 0, 0;\
              -1, 0, 1, 0, 0, 0, 1,-1, 0, 0;\
               0, 0, 0,-1,-1,-1, 0, 0, 0, 0;\
               0, 0, 0, 2, 0, 1,-1, 0, 0, 1;\
               0, 0, 0, 0, 0, 0, 0, 0,-1,-1;\
               0, 0, 0, 0, 0,-1, 1,-1, 2,-1;\
               1,-1, 0, 0, 0, 0, 0, 0, 0, 0;\
               0, 0, 0, 0, 0, 0, 0, 1, 0, 0;\
               0, 0, 0, 1, 1, 1, 0, 0, 0, 0;\
               0, 0, 0, 0, 1, 0, 0, 0, 0, 0')
               


species = ['O3','NO','NO2','HCHO',
           'HO2','HO2H','OH','O',
           'HNO3', 'CO', 'H2']
reactions = ['R' + str(i) for i in range(1,A.shape[1]+1)]
               
# Get indices (row, col) and values for the edgelist               
rows, cols = A.nonzero()
vals = np.array(A[rows,cols])[0]
elist = [(species[rows[i]], reactions[cols[i]], -vals[i]) for i in range(0,len(rows))]

for i in range(0,len(rows)):
    if elist[i][2]<0:
        elist[i] = (elist[i][1],elist[i][0],-elist[i][2])

# Get unique identifiers for the two nodal sets
spc_nodes = np.unique(rows)
rxn_nodes = np.unique(cols + np.max(rows) + 1)

# Construct graph from edgelist
B = nx.DiGraph()
B.add_nodes_from(species, bipartite=0)
B.add_nodes_from(reactions, bipartite=1)
B.add_weighted_edges_from(elist)


# nx.draw(B,pos=nx.random_layout(B,seed=5))
# nx.draw_kamada_kawai(B,with_labels=True)

# Project to unipartite species graph
G = bipartite.projected_graph(B, species,multigraph=True)
# nx.draw_kamada_kawai(G,with_labels=True)

elistG = list(G.edges())
edgesrG = list(G.edges) # with the reactions in index [2] included
unique = [elistG[0]]
duplicates = []
for i in range(1,len(elistG)):
    if elistG[i] not in unique:
        unique.append(elistG[i])
    else:
        duplicates.append(elistG[i])
            

        
# Can we find  cycles?
# note there are other algorithms in this library
B_cycles = sorted(nx.simple_cycles(B))
G_cycles = sorted(nx.simple_cycles(G))
# Note the ['NO2', 'NO', 'OH', 'HO2'] cycle is not a simple cycle in B
# as it passes through R7 twice in one cycle

# Instantaneous concentrations [ppm] and rxn rates r [ppm/min]
c = [0.011074633435759564, 0.0036884511412045525, 0.0019041132919834663,
     0.0020982031519012033, 3.0368122395333126e-6, 0.001343693090096041, 
     2.371154775847784e-6,2.280694342039392e-10, 0.0036848964071020743,
     0.011964393828162825, 0.0024447112659345817]
r = [0.0009520566459917331, 0.0009520566459917331, 0.0008588897336756119,
     3.147304727851805e-5, 4.616046934182647e-5, 6.289127033046298e-5,
     0.00013317892338188386, 6.399884921596744e-5, 4.0310792702881226e-7,
     7.09501968951105e-6]

# rk can be used as a mass rate law double check, does r = rk*[reactants]
rk = [0.5, 21.64992249920054, 25.35720808830782,
      0.015, 0.022, 13546.069967078734, 12468.921210108543,
      15146.818480299688, 0.0003, 2426.0437477352893]

c = np.array(c)*1e3  # probably better to convert to molec/cm^3 from ppb
r = np.array(r)*1e3 

# r = [int(np.array(r[i])*1e6) for i in range(0,len(r))]


seed = 10 # 13640  # Seed random number generators for reproducibility
# G = nx.random_k_out_graph(11, 3, 0.5, seed=seed)
G = bipartite.projected_graph(B, species,multigraph=True)
pos = nx.random_layout(G, seed=seed)
# Tweak the positions, lol maybe just hard code them next time
pos["H2"] += (-0.80, -.02)
pos["HO2H"] += (-0.7, -0.5)
pos["HO2"] += (-0.2, -0.01)
pos["CO"] += (-0.95, -0.04)
pos["OH"] += (0.2, -0.1)
pos["HNO3"] += (-0.35, 0.00)
pos["NO"] += (-0.15, -0.2)
pos["NO2"] += (0.15, 0.02)
pos["HCHO"] += (-0.4, -0.2)



# r_new = []
# for reactant in species: 
#     for i in range(0,len(edgesrG)-1):
#         if reactant == edgesrG[i][0]: 
#             r_new.append((-1)*r[int(edgesrG[i][2][1])-1])
#         elif reactant == edgesrG[i][1]: # for a product
#             r_new.append(r[int(edgesrG[i][2][1])-1])
            
    
edgesrG = list(G.edges)
r_edges = np.array([ r[int(edge[2][1])-1] for edge in G.edges])
simple_edges = list(nx.Graph(G).edges())
r_simple = np.zeros(len(simple_edges))
for i in range(0,len(r_simple)):
    for edge in G.edges:
        if set(edge[0:2]) == set(simple_edges[i]):
            if edge[0:2] == simple_edges[i]:
                r_simple[i] += r[int(edge[2][1])-1] 
            else:
                r_simple[i] -= r[int(edge[2][1])-1] 
    if r_simple[i] < 0:
        simple_edges[i] = (simple_edges[i][1],simple_edges[i][0])
        r_simple[i] = -r_simple[i]       
            
        


# r_edges = np.array([ r[int(edge[2][1])-1] for edge in G.edges])

node_sizes = 30*np.log10(c*1e9) # just so every decade is an order of magnitude
M = len(r_simple) # previously G.number_of_edges() prior to net flow
edge_colors = r_simple #used to be r_edges  # r/(np.max(r)+1e-2)# range(2, 35)# 
edge_alphas =  1+np.log10(r_simple)/5
cmap = plt.cm.plasma

edgelistnew = []
for i in range(0, len(edgesrG)-1):
    edgelistnew.append(list(edgesrG[i]))
    edgelistnew[i].append('rad=0.1')
    

nodes = nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color="lightblue")
nx.draw_networkx_labels(G, pos,font_size=6)
# # edges = nx.draw_networkx_edges(
# #     G,
# #     pos,
# #     node_size=node_sizes,
# #     arrowstyle="->",
# #     arrowsize=10,
# #     edge_color=edge_colors,
# #     edge_vmin = np.min(edge_colors),
# #     edge_vmax = np.max(edge_colors),
# #     edge_cmap=cmap,
# #     width=2,
# #     # connectionstyle='arc3, rad=rrr'.replace('rrr',str(0.3*e[2])),
# #     # connectionstyle='arc3, rad=rrr'.replace('rrr',str(0.3*(r_edges[i]) for i in range(0, len(r_edges)-1)))
# #     )

# for i in range(M):
#     edges[i].set_alpha(edge_alphas[i])



# pc = mpl.collections.PatchCollection(edges, cmap, norm=mpl.colors.LogNorm())
# pc.set_array(edge_colors)

# ax = plt.gca()
# ax.set_axis_off()
# plt.colorbar(pc, ax=ax)
# plt.show()


edge_color = edge_colors
# nx.draw(MDG, pos, with_labels=True, edge_color = (1,1,1))
# for e in G.edges:
# edgesdrawn = []
for e, ec in zip(simple_edges, edge_color):
    color = cmap(ec)
    plt.gca().annotate("",
                       xy=pos[e[1]], 
                       xycoords='data',
                       xytext=pos[e[0]], 
                       textcoords='data',
                       arrowprops=dict(arrowstyle="->", color=color, # edge_color=edge_colors, 
                                       # edge_vmin = np.min(edge_colors), edge_vmax = np.max(edge_colors), edge_cmap=cmap,
                                       shrinkA=10, shrinkB=10,
                                       patchA=None, patchB=None,
                                       connectionstyle="arc3,rad=rrr".replace('rrr',str(rd.random()*0.5+0.1)))
                      )
    # edgesdrawn.append(e)

norm = mpl.colors.Normalize(vmin=np.min(edge_colors), vmax=np.max(edge_colors))
# pc = mpl.collections.PatchCollection(edgesdrawn, cmap, norm=mpl.colors.LogNorm())
pc = plt.cm.ScalarMappable(cmap=cmap, norm=mpl.colors.LogNorm())
pc.set_array(edge_colors)

plt.axis('off')
plt.colorbar(pc)
# plt.show()

# edges = []
# for i in range(0,len(r_simple)-1):
#     edges.append(nx.draw_networkx_edges(
#         G,
#         pos,
#         node_size=node_sizes,
#         arrowstyle="->",
#         arrowsize=10,
#         edge_color=edge_colors,
#         edge_vmin = np.min(edge_colors),
#         edge_vmax = np.max(edge_colors),
#         edge_cmap=cmap,
#         width=2,
#         connectionstyle='arc3, rad=rrr'.replace('rrr',str(0.3*e[2])),
#         connectionstyle='arc3, rad=rrr'.replace('rrr',str(0.3*i))
#         ))


# set alpha value for each edge
# for edge in G.edges(data=True):
#     nx.draw_networkx_edges(G, pos, connectionstyle=f'arc3, rad = {edge[2]["rad"]}')


plt.savefig('/Users/emywli/Desktop/Research/MultigraphGradientNetFlow.pdf')




