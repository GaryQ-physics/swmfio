import swmfio

import logging
swmfio.logger.setLevel(logging.INFO)

import os
import glob
files = glob.glob('/Volumes/My Passport for Mac/git-data/dwelling/divB_simple1/GM/3d__*.out')
for file in files:
    if os.path.exists(file[0:-4] + ".vtk"):
        print("Skipping " + file)
        continue
    batsclass = swmfio.read_batsrus(file)
    vtkfile = swmfio.write_vtk(batsclass)

