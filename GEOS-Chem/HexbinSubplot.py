import matplotlib.pyplot as plt
import numpy as np

# Generate random uniform data
x1 = np.random.rand(1000)
y1 = np.random.rand(1000)

# Generate random Gaussian data
x2 = np.random.randn(10000)
y2 = np.random.randn(10000)

# Create a figure and two subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4), sharey=True)

# Plot the first hexbin
ax1.hexbin(x1, y1, gridsize=20, cmap='Blues')
ax1.set_title('Hexbin Plot 1')

# Plot the second hexbin
ax2.hexbin(x2, y2, gridsize=30, cmap='Reds')
ax2.set_title('Hexbin Plot 2')

# Add a colorbar
cbar = fig.colorbar(ax1.collections[0], ax=[ax1], orientation='horizontal')
cbar.set_label('Counts')

# Add another colorbar
cbar = fig.colorbar(ax2.collections[0], ax=[ax2], orientation='horizontal')
cbar.set_label('Counts')    

# Show the plot
plt.show()