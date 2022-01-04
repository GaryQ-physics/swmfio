import numpy as np
from os.path import exists
from urllib.request import urlretrieve
from swmf_file_reader.batsrus_class import get_class_from_native

urlbase = 'http://mag.gmu.edu/git-data/swmf_file_reader/demodata/'
tmpdir = '/tmp/'
filebase = '3d__var_2_e20190902-041000-000'

for ext in ['.tree', '.info', '.out']:
    filename = tmpdir + filebase + ext
    if not exists(filename):
        print("Downloading " + urlbase + filebase + ext)
        urlretrieve(urlbase + filebase + ext, tmpdir + filebase + ext)

# Instantiate
print("Reading " + tmpdir + filebase + ".*")
batsclass = get_class_from_native(tmpdir + filebase)

# Print methods and attributes
print(dir(batsclass))

# Get data on native grid
print( batsclass.data_arr.shape )
# (5896192, 19)

# Interpolate
print( batsclass.interpolate(np.array([1.,1.,1.]), 'rho') )
# 9.5159912109375

# Derived quantities
print( batsclass.get_native_partial_derivatives(123456, 'rho') )
# [0.25239944 0.41480255 0.7658005 ]