import numpy as np
from os.path import exists
from urllib.request import urlretrieve
from swmf_file_reader.batsrus_class import get_class_from_cdf

urlbase = 'http://mag.gmu.edu/git-data/swmf_file_reader/demodata/'
tmpdir = '/tmp/'
filename = '3d__var_2_e20190902-041000-000.out.cdf'

if not exists(tmpdir + filename):
    print("Downloading " + urlbase + filename)
    print("to")
    print(tmpdir + filename)
    urlretrieve(urlbase + filename, tmpdir + filename)

batsclass = get_class_from_cdf(tmpdir + filename)

assert batsclass.data_arr.shape == (5896192, 19)

# Get a 513th value of x, y, z, and rho
var_dict = dict(batsclass.varidx)
rho = batsclass.data_arr[:, var_dict['rho']][513]
x = batsclass.data_arr[:,var_dict['x']][513]
y = batsclass.data_arr[:,var_dict['y']][513]
z = batsclass.data_arr[:,var_dict['z']][513]
print(f'x={x}, y={y}, z={z}; rho={rho}')
# x=-71.25, y=-15.75, z=-7.75; rho=4.228740215301514

# Get interpolated value at a native grid point
rhoi = batsclass.interpolate(np.array([x, y, z]), 'rho')

assert rho == rhoi

print( batsclass.get_native_partial_derivatives(123456, 'rho') )
# [0.25239944 0.41480255 0.7658005 ]