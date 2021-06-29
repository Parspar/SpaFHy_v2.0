# -*- coding: utf-8 -*-
"""
Created on Fri Jun  4 13:47:49 2021

@author: janousu
"""

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import numpy as np
import matplotlib.colors as mcolors
import numpy as np

# reading the results
outputfile = 'C:\SpaFHy_v1_Pallas_2D/results/testcase_input_d4.nc'
results = read_results(outputfile)

xcrop = np.arange(24,162)
ycrop = np.arange(21,245)

x_plain = np.zeros(results['parameters_cmask'][ycrop,xcrop].shape)
y_plain = np.zeros(results['parameters_cmask'][ycrop,xcrop].shape)
z_plain = results['soil_ground_water_level'][20,ycrop,xcrop]
z_plain2 = results['parameters_elevation'][ycrop, xcrop] - results['soil_ground_water_level'][20,ycrop,xcrop]
dem_plain = results['parameters_elevation'][ycrop, xcrop]


'''
x_plain = np.rot90(x_plain, k=1, axes=(0, 1))
y_plain = np.rot90(y_plain, k=1, axes=(0, 1))
z_plain = np.rot90(z_plain, k=1, axes=(0, 1))
dem_plain = np.rot90(dem_plain, k=1, axes=(0, 1))
'''

for i in range(x_plain.shape[0]):
    for j in range(x_plain.shape[1]):
        i_coord = i*1
        j_coord = j*1
        x_plain[i,j] = i_coord
        y_plain[i,j] = j_coord
        
        
#%%

norm = mcolors.TwoSlopeNorm(vmin=np.nanmin(z_plain), vmax = np.nanmax(z_plain), vcenter=0)
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

surf = ax.plot_surface(x_plain, y_plain, z_plain, cmap=cm.coolwarm_r,norm=norm,
                       linewidth=0, antialiased=False)

# Customize the z axis.
ax.set_zlim(-40, 40)
# A StrMethodFormatter is used automatically
ax.zaxis.set_major_formatter('{x:.02f}')

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)
#ax4.axis("off")
ax.view_init(40, 225)


plt.show()


#%%

norm = mcolors.TwoSlopeNorm(vmin=np.nanmin(dem_plain), vmax = np.nanmax(dem_plain), vcenter=np.nanmean(dem_plain))
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})


surf2 = ax.plot_surface(x_plain, y_plain, dem_plain,cmap='terrain',norm=norm,
                       linewidth=0, antialiased=False)
# Customize the z axis.
ax.set_zlim(-5, 600)
# A StrMethodFormatter is used automatically
ax.zaxis.set_major_formatter('{x:.02f}')

# Add a color bar which maps values to colors.
fig.colorbar(surf2, shrink=0.5, aspect=5)

ax.view_init(40, 225)

plt.show()


#%%

# together


norm = mcolors.TwoSlopeNorm(vmin=np.nanmin(z_plain), vmax = np.nanmax(z_plain), vcenter=0)
norm2 = mcolors.TwoSlopeNorm(vmin=np.nanmin(dem_plain), vmax = np.nanmax(dem_plain), vcenter=np.nanmean(dem_plain))

fig = plt.figure(figsize=plt.figaspect(0.5))

#fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

ax = fig.add_subplot(1, 2, 1, projection='3d')

#ax1 = ax[0]
#ax2 = ax[1]

surf = ax.plot_surface(x_plain, y_plain, z_plain, cmap=cm.coolwarm_r,norm=norm,
                       linewidth=0, antialiased=False)

# Customize the z axis.
ax.set_zlim(np.nanmin(z_plain)-50, np.nanmax(z_plain)+50)
#ax.set_zlim(-40, 40)
# A StrMethodFormatter is used automatically
ax.zaxis.set_major_formatter('{x:.02f}')
ax.title.set_text('Groundwater depth [m]')
# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)
#ax4.axis("off")
ax.view_init(30, 210)

ax = fig.add_subplot(1, 2, 2, projection='3d')

surf2 = ax.plot_surface(x_plain, y_plain, dem_plain,cmap='terrain',norm=norm2,
                       linewidth=0, antialiased=False)
ax.title.set_text('Surface elevation [m]')

# Customize the z axis.
ax.set_zlim(np.nanmin(dem_plain)-300, np.nanmax(dem_plain)+300)
#ax.set_zlim(100, 500)
# A StrMethodFormatter is used automatically
ax.zaxis.set_major_formatter('{x:.02f}')

# Add a color bar which maps values to colors.
fig.colorbar(surf2, shrink=0.5, aspect=5)
#ax4.axis("off")
ax.view_init(30, 210)



plt.show()



