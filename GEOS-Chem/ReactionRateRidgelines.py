import numpy as np
import networkx as nx
import pandas as pd
import seaborn as sns
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from ChemicalCase import ChemicalCase
import os

# Important: keep this set_theme call for Emy-style Plots in the KDEs
sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})



# Read the file into a DataFrame
df_species = pd.read_csv('../SPC_NAMES.txt', header=None, delim_whitespace=True)

# find location of select species in the dataframe
# species_list = ["ISOP", "MTPA", "O3", "NO", "HO2"]
species_list = ["O3", "NO", "OH", "HO2", "Cl", "ClO"]
# colors are green, dark green, lime green, grey, brown
species_colors = ['g', 'darkgreen', 'grey', 'brown']
species_loc = np.ones(len(species_list), dtype=int) * -1
for i, species in enumerate(species_list):
    species_loc[i] = df_species[df_species[0] == species].index[0]
species_loc = species_loc.tolist()


# for each location, get the chemical case 
# from the file in "L1 Noon Samples" folder
samples_directory = '../L1 Noon Samples/'

# Add the 16 location names
loc_names = ['Atlantic Ocean', 'Indian Ocean', 'Pacific Ocean', 'Graciosa', 'Kennaook', 
             'El Djouf', 'Amazon', 'Borneo', 'Congo', 'Ozarks', 
             'Beijing', 'Kinshasa', 'Los Angeles', 'Paris', 
             'Utqiagvik', 'McMurdo Station']
# loc_names = ['Los Angeles', "Beijing", "Paris", "Ozarks"]
# all forest locations
# loc_names = ['Amazon', 'Borneo', 'Congo', 'Ozarks']
file_prefix = [name.replace(" ","").replace("Station","").replace("Kennaook","CapeGrim") for name in loc_names]

# find the file that starts with the file prefix
file_names = []
for prefix in file_prefix:
    for file in os.listdir(samples_directory):
        if file.startswith(prefix):
            file_names.append(file)
            break

cases = []
for file in file_names:
    cases.append(ChemicalCase(samples_directory + file))

loc_dict = {
    "Atlantic Ocean": [0.2154516905067, 0.295, 0.5078367867510504],
    "Indian Ocean": [0.2154516905067, 0.295, 0.5078367867510504],
    "Pacific Ocean": [0.2154516905067, 0.295, 0.5078367867510504],
    "Graciosa": [0.2154516905067, 0.295, 0.5078367867510504],
    "Kennaook": [0.2154516905067, 0.295, 0.5078367867510504],
    "El Djouf": [0.78, 0.51, 0.1892045842],
    "Amazon": [0.12071162840208301, 0.4357915995132193, 0.2463679091477368],
    "Borneo": [0.12071162840208301, 0.4357915995132193, 0.2463679091477368],
    "Congo": [0.12071162840208301, 0.4357915995132193, 0.2463679091477368],
    "Ozarks": [0.12071162840208301, 0.4357915995132193, 0.2463679091477368],
    "Beijing": [0.7, 0.7, 0.65],
    "Kinshasa": [0.7, 0.7, 0.65],
    "Los Angeles": [0.7, 0.7, 0.65],
    "Paris": [0.7, 0.7, 0.65],
    "Utqiagvik": [0.5184891135337089, 0.7194769933793438, 0.7518803726887796],
    "McMurdo Station": [0.5632111255041908, 0.758620966612444, 0.7764634182455044]
}
loc_colors = [loc_dict[loc] for loc in loc_names]

# convert cases object to a dataframe, where each row is a reaction rate
# an additional "tiles" column is added to store the location name
log_reaction_rates = []
tiles = []
for i, case in enumerate(cases):
    for rate in case.reaction_rates:
        # skip if rate is <= 1e-10
        if rate <= 1e-25:
            continue
        log_reaction_rates.append(np.log10(rate))
        tiles.append(loc_names[i])

# Create DataFrame
df = pd.DataFrame(dict(reaction_rate=log_reaction_rates, tiles=tiles))

# Color palette, each category of location having a different color
pal = sns.cubehelix_palette(16, rot=-.25, light=.7)
pal[15] = pal[0]
pal[14] = pal[1]
pal[0] = pal[1] = pal[2] = pal[3] = pal[4] = [0.2154516905067, 0.295, 0.5078367867510504]
pal[5] = [0.78, 0.51, 0.1892045842]
pal[6] = pal[7] = pal[8] = pal[9] = [0.12071162840208301, 0.4357915995132193, 0.2463679091477368]
pal[10] = pal[11] = pal[12] = pal[13] = [0.7, 0.7, 0.65]

tiles = sns.FacetGrid(df, row="tiles", hue="tiles", aspect=15, height=.5, palette=pal)

# Mapping the densities
tiles.map(sns.kdeplot, "reaction_rate",
      bw_adjust=.5, clip_on=False,
      fill=True, alpha=1, linewidth=1.5)
tiles.map(sns.kdeplot, "reaction_rate", clip_on=False, color="w", lw=2, bw_adjust=.5)

# Passing color=None to refline() uses the hue mapping
tiles.refline(y=0, linewidth=2, linestyle="-", color=None, clip_on=False)

def label(x, color, label):
    ax = plt.gca()
    ax.text(0, .2, label, fontweight="bold", color=color,
            ha="left", va="center", transform=ax.transAxes)
    ax.set_xlim(-25, 10)

tiles.map(label, "reaction_rate")

# Set the subplots to overlap
tiles.figure.subplots_adjust(hspace=-0.15)

# Remove axes details that don't play well with overlap
tiles.set_titles("")

ticks = plt.xticks()
xtick_string = [rf'$10^{str({i})}$' for i in range(-25, 11, 5)]

tiles.set(xticks=range(-25, 11, 5))
tiles.set_xticklabels(xtick_string)

tiles.set(xlabel='Reaction Rate [molec/cm^3/s]')
tiles.set(yticks=[], ylabel="")

tiles.despine(bottom=True, left=True)

plt.show()