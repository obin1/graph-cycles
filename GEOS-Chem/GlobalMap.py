import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.patheffects as PathEffects
import cartopy.crs as ccrs
import os
import ChemicalCase
import cartopy.feature as cfeature



# Create a figure and axis
fig = plt.figure(figsize=(10, 6))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())


# Set the extent to the whole world
ax.set_global()

# Add satellite imagery
ax.stock_img()


# Add the 16 location names
loc_names = ['Atlantic Ocean', 'Indian Ocean', 'Pacific Ocean', 'Graciosa', 'Kennaook', 
             'El Djouf', 'Amazon', 'Borneo', 'Congo', 'Ozarks', 
             'Beijing', 'Kinshasa', 'Los Angeles', 'Paris', 
             'Utqiagvik', 'McMurdo Station']

loc_name_lat_offsets = [-3.5]*len(loc_names)
loc_name_lon_offsets = [-5]*len(loc_names)
# And hard code a few locations to look nice
# Graciosa needs to be moved up and to the right
loc_name_lat_offsets[3] = 3
loc_name_lon_offsets[3] = 5
# Atlantic Ocean needs to be moved down and to the right
loc_name_lat_offsets[0] = -10
loc_name_lon_offsets[0] = 5
# Ozarks needs to be moved up and to the right
loc_name_lat_offsets[9] = 3
loc_name_lon_offsets[9] = 5
# Utqiagvik needs to reversed so that it is on the right side of the point
loc_name_lon_offsets[14] = 35
# Congo needs to be reversed so that it is on the right side of the point
loc_name_lon_offsets[8] = 26
# El Djouf needs to be reversed so that it is on the right side of the point
loc_name_lon_offsets[5] = 30



loc_colors = [[0.2154516905067, 0.295, 0.5078367867510504],  # Atlantic Ocean
              [0.2154516905067, 0.295, 0.5078367867510504],  # Indian Ocean
              [0.2154516905067, 0.295, 0.5078367867510504],  # Pacific Ocean
              [0.2154516905067, 0.295, 0.5078367867510504],  # Graciosa
              [0.2154516905067, 0.295, 0.5078367867510504],  # Kennaook
              [0.78, 0.51, 0.1892045842],  # El Djouf
              [0.12071162840208301, 0.4357915995132193, 0.2463679091477368],  # Amazon
              [0.12071162840208301, 0.4357915995132193, 0.2463679091477368],  # Borneo
              [0.12071162840208301, 0.4357915995132193, 0.2463679091477368],  # Congo
              [0.12071162840208301, 0.4357915995132193, 0.2463679091477368],  # Ozarks
              [0.7, 0.7, 0.65],  # Beijing
              [0.7, 0.7, 0.65],  # Kinshasa
              [0.7, 0.7, 0.65],  # Los Angeles
              [0.7, 0.7, 0.65],  # Paris
              [0.5184891135337089, 0.7194769933793438, 0.7518803726887796],  # Utqiagvik
              [0.5632111255041908, 0.758620966612444, 0.7764634182455044]]  # McMurdo Station
file_names_notwilight = [name.replace(" ","").replace("Station","").replace("Kennaook","CapeGrim")+"_L1_20180701_2100.txt" for name in loc_names]
file_names_twilight = [name.replace(" ","").replace("Station","").replace("Kennaook","CapeGrim")+"Twilight_L1_20180701_2100.txt" for name in loc_names]
file_names_combined = [None]*(len(file_names_notwilight)+len(file_names_twilight))
file_names_combined[::2] = file_names_notwilight
file_names_combined[1::2] = file_names_twilight

# Directory to search
samples_directory = 'samples/'

# Initialize an empty list for file names
file_names = []

# Make the actual list of file names
for file_name in file_names_combined:
    if os.path.isfile(os.path.join(samples_directory, file_name)):
        file_names.append(file_name)

# Initialize a list of ChemicalCase objects
cases = []
for sample in file_names:
    cases.append(ChemicalCase.ChemicalCase(samples_directory + sample))

# Plot the locations
for i in range(len(cases)):
    ax.scatter(cases[i].longitude, cases[i].latitude, marker='*', 
                   color=loc_colors[i], s=100, edgecolors='black', linewidth=0.4,
                   transform=ccrs.PlateCarree())
    # add a text box against a black background
    txt = ax.text(cases[i].longitude+loc_name_lon_offsets[i], cases[i].latitude+loc_name_lat_offsets[i], 
            loc_names[i], color = loc_colors[i], transform=ccrs.PlateCarree(),
            horizontalalignment='right', verticalalignment='bottom',backgroundcolor='none')
    txt.set_path_effects([PathEffects.withStroke(linewidth=0.7, foreground='black')])
plt.draw()  

plt.savefig('GlobalMap.png', dpi=300)
