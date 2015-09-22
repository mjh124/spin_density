# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 15:56:11 2015

@author: michael
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import DataFileParser

# Script to plot spin density cross section of an arbitrary plane (defined by 3 points) from a cube file

MyData = '../Desktop/Co13PH3Cl6_half_spindens.cube'
param, mass, x, y, z, dens = DataFileParser.ParseCube(MyData)

scale = 0.377945 # Width of each voxel (bohr)
Natoms = param[0]
xgrid = param[1]
ygrid = param[2]
zgrid = param[3]

# Find (or define) the 3 atom positions and their coordinates used to define the plane

# Define the three points
p1 = np.array([x[0]/scale,y[0]/scale,z[0]/scale])
p2 = np.array([x[17]/scale,y[17]/scale,z[17]/scale])
p3 = np.array([x[13]/scale,y[13]/scale,z[13]/scale])
# Define vectors in the plane
p1p2 = p2 - p1 # vector from center Co to terminal Cl
p1p2_len = np.linalg.norm(p1p2) # length of vector from Co to terminal Cl
p1p3 = p3 - p1 # vector from center Co to other Cl
# Determine normal vector to plane
norm = np.cross(p1p2, p1p3)
norm_len = np.linalg.norm(norm)
print norm
a, b, c = norm
# Slove for d
d = np.dot(norm, p3)
# Write equation of the plane
print('The equation of the plane is {0}x + {1}y + {2}z = {3}'.format(a, b, c, d))

# Plot the plane
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#x = np.linspace(0, xgrid, xgrid-1)
#y = np.linspace(0, ygrid, ygrid-1)
#X, Y = np.meshgrid(x, y)
#Z = (d - a*X - b*Y)/c
## Plot the meshgrid
#ax.plot(X.flatten(),
#        Y.flatten(),
#        Z.flatten(), 'bo ', markersize=1)
## Plot the original points
#ax.plot(*zip(p1, p2, p3), color='r', linestyle=' ', marker='o')
## Adjust the view so we can see the point/plane alignment
##ax.view_init(0,22)
#plt.tight_layout()
##plt.savefig('../Desktop/plane.png')
#plt.show()

# Constants to be used as conversion factors
gp_to_ang = 0.20 # Convert grip point indices to angstroms
rad_to_deg = 180/np.pi # Convert radians to degrees

# Setup binning system
Nbin_dist = 250.0
bin_dist_max = 25.0
bin_dist_size = bin_dist_max / Nbin_dist

Nbin_ang = 360.0
bin_ang_max = 360.0
bin_ang_size = bin_ang_max / Nbin_ang

# Calculate the projection of a vector from center Co to any point onto the plane
#radial = []
#angular = []
Ngrid = xgrid*ygrid*zgrid
map = np.zeros((Nbin_ang, Nbin_dist))
for i in range(Ngrid):
    
    rho = dens[i]
    
    xc = (i // (ygrid*zgrid))
    yc = ((i // zgrid) % ygrid)
    zc = (i % zgrid)
    
    p = np.array([xc, yc, zc]) - p1 # Vector from center Co to point p

    # Project the point p onto the plane defined earlier    
    dot = np.dot(p, norm)
    proj_n = (dot/(np.linalg.norm(norm))**2) * norm # Calculate the projection of p onto the normal of the plane
    norm_comp = np.linalg.norm(proj_n) # Calculate length of projection onto the normal vector for filter
    projection = p - proj_n # Take difference of p and proj_n to get component in plane

    if norm_comp >= 2.0:
        continue
    else:

        # Measure the length of the projection vector
        proj_len = np.linalg.norm(projection) # Units of grid points
        proj_len_ang = proj_len * gp_to_ang
        col = int(proj_len_ang // bin_dist_size)
        #    radial.append(proj_len_ang)

        # Calculate the cross product of the projection and the Co-Cl vector for comparison with norm
        cross = np.cross(projection, p1p2)
        cross_len = np.linalg.norm(cross)
        tmp = np.dot(cross, norm) / (cross_len*norm_len)
        cross_norm_ang = np.arccos(tmp)      

        # Measure the angle between the projection and the Co-Cl vector
        tmp1 = np.dot(projection, p1p2) / (proj_len*p1p2_len)
        ang = np.arccos(tmp1)

        if cross_norm_ang > 3.13:
            angle = (2*np.pi - ang)
            angle_deg = angle * rad_to_deg
            row = int(angle_deg // bin_ang_size)

        else:
            angle_deg = ang * rad_to_deg
            row = (angle_deg // bin_ang_size)
                
        map[row, col] = map[row, col] + rho

print map.shape
#for i in range(35):
#    for j in range(360):
#        print i, j, map[i][j]
#print min(radial), max(radial)
#print min(angular), max(angular)

# Plot the radial heatmap

azimuths = np.radians(np.linspace(0, bin_ang_max, Nbin_ang, endpoint=True))
#azimuths = np.radians(np.linspace(0, 180, 180, endpoint=False))
zeniths = np.linspace(0, bin_dist_max, Nbin_dist, endpoint=True)

r, theta = np.meshgrid(zeniths, azimuths)

#x = r*np.cos(theta)
#y = r*np.sin(theta)

fig = plt.figure()
cmap = plt.get_cmap("hot")
ax = fig.add_subplot(111, polar = True)
CS = ax.pcolor(theta, r, map, cmap=cmap)
ax.set_rmax(6)
#image = plt.imshow(cmap=cmap, interpolation=gaussian)

#data = griddata()

CB = ax.colorbar(CS, shrink=0.8, extend='both')
plt.show()
#plt.savefig('../Desktop/Co13PH3Cl6_half_SDCS_360.png', bbox_inches='tight')
#plt.clf()

#CS = plt.contour(x, y, map)
#plt.show()