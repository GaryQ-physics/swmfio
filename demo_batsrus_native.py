from timeit import default_timer as timer
from datetime import timedelta

import swmfio as swmfio

if True:
    # Read data in .out file into numpy array with columns x, y, z, ...
    # Column names are given in vars
    # No tree information is returned.
    #url = 'http://mag.gmu.edu/git-data/dwelling/divB_simple1/GM/3d__mhd_4_e20100320-013000-000'
    url = 'http://mag.gmu.edu/git-data/swmfio/3d__var_2_e20190902-041000-000'
    file = swmfio.dlfile(url, progress=True)

    from swmfio.read_batsrus import read_data
    data, vars, meta = read_data(file)
    print("data")
    print(data)
    print("vars")
    print(vars)
    print("meta")
    print(meta)

if True:
    urlbase = 'http://mag.gmu.edu/git-data/swmfio/3d__var_2_e20190902-041000-000'
    filebase = swmfio.dlfile(urlbase, progress=True)

    start = timer()
    print("Reading and creating native grid interpolator for " + filebase + ".{tree, info, out}")
    batsclass = swmfio.read_batsrus(filebase)
    end = timer()
    print("Time: {}".format(timedelta(seconds=end-start)))
    # Time: 0:00:01.856912

    assert batsclass.data_arr.shape == (5896192, 20)

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
    rhoi = batsclass.interpolate([x, y, z], 'rho')
    assert rho == rhoi
    end = timer()
    print("Interpolation time: {}".format(timedelta(seconds=end-start)))
    # Interpolation time: 0:00:00.004919
    assert rho == rhoi

    # Compute a derived quantity
    print( batsclass.get_native_partial_derivatives(123456, 'rho') )
    # [0.25239944 0.41480255 0.7658005 ]

    if False:
        # Demonstrate/check interpolation without swmfio. This takes a very long time
        # (> 10 minutes if it does not throw error)
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
        # Creation time: 0:11:47.466515

        start = timer()
        print("Interpolating using LinearNDInterpolator")
        rhoi = rho_interpolator([x, y, z])
        end = timer()
        print("Interpolation time: {}".format(timedelta(seconds=end-start)))
        # Interpolation time: 0:00:58.665182
        assert rho == rhoi
