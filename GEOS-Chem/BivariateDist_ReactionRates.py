import numpy as np
import networkx as nx
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from ChemicalCase import ChemicalCase
import os

# obtaining species and converting to list from SPC_NAMES.txt
with open('SPC_NAMES.txt', "rt") as my_file:
    species = my_file.read().split()

# for each location, get the chemical case 
# from the file in "L1 Noon Samples" folder
# samples_directory = 'L1 Noon Samples/'
samples_directory = 'samples/'

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

# make a dictionary of location types
loc_type_dict = {
    "AtlanticOcean": "ocean",
    "IndianOcean": "ocean",
    "PacificOcean": "ocean",
    "Graciosa": "ocean",
    "CapeGrim": "ocean",
    "ElDjouf": "desert",
    "Amazon": "biogenic",
    "Borneo": "biogenic",
    "Congo": "biogenic",
    "Ozarks": "biogenic",
    "Beijing": "urban",
    "Kinshasa": "urban",
    "LosAngeles": "urban",
    "Paris": "urban",
    "Utqiagvik": "arctic",
    "McMurdo": "antarctic"
}

loc_colors_dict = {
    "ocean": [0.2154516905067, 0.295, 0.5078367867510504],
    "desert": [0.78, 0.51, 0.1892045842],
    "biogenic": [0.12071162840208301, 0.4357915995132193, 0.2463679091477368],
    "urban": [0.7, 0.7, 0.65],
    "arctic": [0.5184891135337089, 0.7194769933793438, 0.7518803726887796],
    "antarctic": [0.5632111255041908, 0.758620966612444, 0.7764634182455044]
}
loc_colors = [loc_colors_dict[loctype] for loctype in list(loc_colors_dict.keys())]


cases = []
for file in file_names:
    cases.append(ChemicalCase(samples_directory + file))

# make a dataframe of all the reaction rates, incuding other things like location etc.

# get locations of NO and NO2 in the concentrations list
spc_names = pd.read_csv('../SPC_NAMES.txt', header=None, delim_whitespace=True)
NO_loc = spc_names[spc_names[0] == 'NO'].index[0]
NO2_loc = spc_names[spc_names[0] == 'NO2'].index[0]


# Initialize a list to store the data
data = []
total_reaction_rate = []
nox_vector = []

for i, case in enumerate(cases):
    # only get surface level 1
    if case.level != 1:
        continue
    # only get daylight hours
    if case.cos_sza < 0.6:
        continue

    # get the NOx concentration
    NOx_concentration = case.concentrations[NO_loc] + case.concentrations[NO2_loc]
    nox_vector.append(NOx_concentration)

    # get the reaction rates
    reaction_rates = case.reaction_rates
    total_reaction_rate.append(sum(reaction_rates))

    # put this into data
    for j, rate in enumerate(reaction_rates):
        data.append(["R"+str(j+1), rate, case.location, loc_type_dict[case.location], NOx_concentration])

# make a dataframe
df = pd.DataFrame(data, columns=["Reaction", "Reaction Rate", "Location", "Location Type", "NOx Concentration"])

# set all zero reaction rates to 1e-50
df.loc[df['Reaction Rate'] == 0, 'Reaction Rate'] = 1e-50

# Plot total reaction rate against NOx concentration
plt.figure(figsize=(10, 6))
# matplotlib scatter plot
plt.scatter(nox_vector, 
            total_reaction_rate, 
            alpha=0.5)
# make log scale in both axes
plt.yscale('log')
plt.xscale('log')
plt.xlabel('NOx Concentration [molec/cm^3]')
plt.ylabel('Total Reaction Rate [molec/cm^3/s]')



# # make a seaborn bivariate kde plot
# plt.figure(figsize=(10, 6))
# sns.kdeplot(
#     data=df,
#     x='NOx Concentration',
#     y='Reaction Rate',
#     # hue='Location Type',
#     fill=True,
#     palette=loc_colors,
#     levels=20,
#     thresh=0.05,
#     alpha=0.5,
#     log_scale=(True, True)
#     # set the x bounds
# )

# find minimum nonzero reaction rate
min_rate = df[df['Reaction Rate'] > 0]['Reaction Rate'].min()