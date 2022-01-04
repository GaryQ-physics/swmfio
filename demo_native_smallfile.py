import numpy as np
from urllib.request import urlretrieve
from swmf_file_reader.batsrus_class import get_class_from_native

urlbase = 'http://mag.gmu.edu/git-data/swmf_file_reader/demodata/'
tmpdir = '/tmp/'
filebase = '3d__mhd_1_t00000000_n00000000'

for ext in ['.tree', '.info', '.out']:
    print("Downloading " + urlbase + filebase + ext)
    urlretrieve(urlbase + filebase + ext, tmpdir + filebase + ext)

# Instantiate
print("Reading " + tmpdir + filebase + ".*")
batsclass = get_class_from_native(tmpdir + filebase)

# Print methods and attributes
print(dir(batsclass))

# Get data on native grid
print( batsclass.data_arr.shape )

# Interpolate
print( batsclass.interpolate(np.array([1.,1.,1.]), 'rho') )

# Derived quantities
print( batsclass.get_native_partial_derivatives(123456, 'rho') )
