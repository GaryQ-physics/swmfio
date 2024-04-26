
###################################################################
###################################################################
#
# Main routine that demos the batsrus interpolation using CCMC
# OCTREE library
#
###################################################################
###################################################################

import swmfio as swmfio
import numpy as np
import random

# Read data in .out file into numpy array with columns x, y, z, ...
# Column names are given in vars
# No tree information is returned.
url = 'http://mag.gmu.edu/git-data/swmfio/3d__var_2_e20190902-041000-000'
file = swmfio.dlfile(url, progress=True)
batsclass = swmfio.read_batsrus(file)
var_dict = dict(batsclass.varidx)

# Demonstrate interpolation using routine borrowed from Kamodo
print( "Creating CCMC AMR interpolator")

# Create an instance of the batsrus interpolator
binterp = swmfio.batsrus_interpolator(batsclass) 

# Register the variables that we will interpolate
binterp.register_variable( 'bx' )
binterp.register_variable( 'by' )
binterp.register_variable( 'bz' )

# We will perform interpolation at some random points in batsrus grid
# Get grid x,y,z
xg = batsclass.data_arr[:,var_dict['x']]
yg = batsclass.data_arr[:,var_dict['y']]
zg = batsclass.data_arr[:,var_dict['z']]

# Get limits of grid
xmin = xg.min()
xmax = xg.max()
ymin = yg.min()
ymax = yg.max()
zmin = zg.min()
zmax = zg.max()

# Create storage for the results
CNT=10
x = np.zeros(CNT)
y = np.zeros(CNT)
z = np.zeros(CNT)
bx = np.zeros(CNT)
by = np.zeros(CNT)
bz = np.zeros(CNT)
bx2 = np.zeros(CNT)

# Pick CNT random points and interpolate the B field at each point
print( f"\nUse interpolator to find Bx, By, Bz at {CNT} random points:\n")
for i in range(CNT):
    x[i] = xmin + random.random() * (xmax-xmin)        
    y[i] = ymin + random.random() * (ymax-ymin)        
    z[i] = zmin + random.random() * (zmax-zmin)
    
    bx[i] = binterp.interp([x[i], y[i], z[i]], 'bx')[0]
    by[i] = binterp.interp([x[i], y[i], z[i]], 'by')[0]
    bz[i] = binterp.interp([x[i], y[i], z[i]], 'bz')[0]

    print(f'{i}: Bx: {bx[i]}, By: {by[i]}, Bz: {bz[i]}')
    
# Demo multiple interpolations in one call
print( f"\nUse interpolator to find Bx at the same {CNT} points, but in one call.")
print( "Bx and Bx2 should be equal:\n")
bx2 = binterp.interp(list(zip(x,y,z)), 'bx')
for i in range(CNT):
    print(f'{i}: Bx: {bx[i]}, Bx2: {bx2[i]}')
    
