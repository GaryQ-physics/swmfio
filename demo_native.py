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
        print("Downloading " + urlbase + filename)
        print("to")
        print(tmpdir + filename)
        urlretrieve(urlbase + filename, tmpdir + filename)

# Instantiate
print("Reading " + tmpdir + filebase + ".*")
batsclass = get_class_from_native(tmpdir + filebase)

# Print methods and attributes
print(dir(batsclass))
# [...,'block2node', 'block_amr_levels', 'block_child_count',
#      'block_child_ids', 'block_parent_id', 'block_x_max',
#      'block_x_min', 'block_y_max', 'block_y_min', 'block_z_max',
#      'block_z_min', 'data_arr', 'find_tree_node',
#      'get_native_partial_derivatives', 'interpolate',
#      'nDim', 'nI', 'nJ', 'nK', 'node2block', 'rootnode',
#      'varidx', 'xGlobalMax', 'xGlobalMin', 'yGlobalMax',
#      'yGlobalMin', 'zGlobalMax', 'zGlobalMin']

assert batsclass.data_arr.shape == (5896192, 19)

# Get a 513th value of x, y, z, and rho
var_dict = dict(batsclass.varidx)
rho = batsclass.data_arr[:, var_dict['rho']][513]
x = batsclass.data_arr[:,var_dict['x']][513]
y = batsclass.data_arr[:,var_dict['y']][513]
z = batsclass.data_arr[:,var_dict['z']][513]
print(x, y, z, rho)
# -71.25 -15.75 -7.75 4.22874

# Get interpolated value at a native grid point
rhoi = batsclass.interpolate(np.array([x, y, z]), 'rho')

assert rho == rhoi

# Compute a derived quantity
print( batsclass.get_native_partial_derivatives(123456, 'rho') )
# [0.25239944 0.41480255 0.7658005 ]