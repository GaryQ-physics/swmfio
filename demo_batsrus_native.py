import numpy as np
from os.path import exists
from urllib.request import urlretrieve

from timeit import default_timer as timer
from datetime import timedelta

import swmfio as swmfio

urlbase = 'http://mag.gmu.edu/git-data/swmfio/'
tmpdir = '/tmp/'
filebase = '3d__var_2_e20190902-041000-000'

for ext in ['.tree', '.info', '.out']:
    filename = filebase + ext
    if not exists(tmpdir + filename):
        print("Downloading " + urlbase + filename)
        print("to")
        print(tmpdir + filename)
        urlretrieve(urlbase + filename, tmpdir + filename)

start = timer()
print("Reading and creating native grid interpolator" + tmpdir + filebase + ".{tree, info, out}")
batsclass = swmfio.read_batsrus(tmpdir + filebase)
end = timer()
print("Time: {}".format(timedelta(seconds=end-start)))
# Time: 0:00:01.856912

assert batsclass.data_arr.shape == (5896192, 19)

# Get a 513th value of x, y, z, and rho
var_dict = dict(batsclass.varidx)
rho = batsclass.data_arr[:, var_dict['rho']][513]
x = batsclass.data_arr[:,var_dict['x']][513]
y = batsclass.data_arr[:,var_dict['y']][513]
z = batsclass.data_arr[:,var_dict['z']][513]
print(x, y, z, rho)
# -71.25 -15.75 -7.75 4.22874

start = timer()
print("Interpolating using batsclass")
rhoi = batsclass.interpolate(np.array([x, y, z]), 'rho')
assert rho == rhoi
end = timer()
print("Interpolation time: {}".format(timedelta(seconds=end-start)))
# Interpolation time: 0:00:00.004919
assert rho == rhoi

if True:
    from scipy.interpolate import LinearNDInterpolator
    start = timer()
    xg = batsclass.data_arr[:,var_dict['x']]
    yg = batsclass.data_arr[:,var_dict['y']]
    zg = batsclass.data_arr[:,var_dict['z']]
    rhog = batsclass.data_arr[:, var_dict['rho']]

    print("Creating LinearNDInterpolator for rho")
    start = timer()
    rho_interpolator = LinearNDInterpolator(list(zip(xg, yg, zg)), rhog)
    end = timer()
    print("Creation time: {}".format(timedelta(seconds=end-start)))

    start = timer()
    print("Interpolating using LinearNDInterpolator")
    rhoi = rho_interpolator(np.array([x, y, z]))
    end = timer()
    print("Interpolation time: {}".format(timedelta(seconds=end-start)))
    assert rho == rhoi

# Compute a derived quantity
print( batsclass.get_native_partial_derivatives(123456, 'rho') )
# [0.25239944 0.41480255 0.7658005 ]