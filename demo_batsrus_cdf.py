import numpy as np
from os.path import exists
from urllib.request import urlretrieve

from timeit import default_timer as timer
from datetime import timedelta

import swmfio as swmfio

urlbase = 'http://mag.gmu.edu/git-data/swmfio/'
tmpdir = '/tmp/'
filename = '3d__var_2_e20190902-041000-000.out.cdf'

if not exists(tmpdir + filename):
    print("Downloading " + urlbase + filename)
    print("to")
    print(tmpdir + filename)
    urlretrieve(urlbase + filename, tmpdir + filename)

start = timer()
batsclass = swmfio.read_batsrus(tmpdir + filename)
end = timer()
print("Read time: {}".format(timedelta(seconds=end-start)))

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
start = timer()
rhoi = batsclass.interpolate(np.array([x, y, z]), 'rho')
assert rho == rhoi
end = timer()
print("Interpolation time: {}".format(timedelta(seconds=end-start)))

# Get interpolated value at a non-native grid point
start = timer()
rhoi = batsclass.interpolate(np.array([x+0.01, y+0.01, z+0.01]), 'rho')
end = timer()
print("Interpolation time: {}".format(timedelta(seconds=end-start)))

print( batsclass.get_native_partial_derivatives(123456, 'rho') )
# [0.25239944 0.41480255 0.7658005 ]