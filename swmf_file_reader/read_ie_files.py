import numpy as np
import re
from numba import types
from numba.typed import Dict


def grep_dash_o(RE, lines):
    ret = ''
    for line in lines:
        findall = re.findall(RE, line)
        if findall != []:
            ret = ret + '\n'.join(findall) + '\n'
    return ret

def read_iono_tec(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()

    var_string = grep_dash_o(r'"\S* \[\S*\]"',lines[1:15])

    varidx = Dict.empty(
        key_type=types.unicode_type,
        value_type=types.int32,
        )
    units = {}
    iVar = np.int32(0)
    for var in var_string.split('\n'):
        if var == '': continue
        split = var[1:-1].split(' ')
        varidx[split[0]] = iVar
        units[split[0]] = split[1][1:-1]
        iVar += np.int32(1)

    assert(len(varidx.keys()) == 22)

    # note, the header has 27 quantities, but last 5 of them have extra white space and are not 
    # picked up by the above grep
    # for even i, the lines will end up to be of size 5 and all zero, so we dont care about them.
    # The first 22 variable names (which are the ones matched by the regex) are then assumed to be the 22 data point on the lines for odd i.

    # note in the SWMF source code, it would appear  Phi <==> Psi  , they are two different names for same thing...
    # likely done at some point to avoid confusion with PHI (the electric potential)
    # this could pose a problem in the future if we wish to convert all keys to lower case...

    north = []
    for i in range(17, 17 + 32942):
        if i%2 == 0:
            assert(np.all( 0 == np.fromstring(lines[i], dtype=np.float32, sep=' ') ))
        elif i%2 == 1:
            north.append(np.fromstring(lines[i], dtype=np.float32, sep=' '))
        else:
            assert(False)

    south = []
    for i in range(17 + 32942 + 2, len(lines)):
        if i%2 == 0:
            assert(np.all( 0 == np.fromstring(lines[i], dtype=np.float32, sep=' ') ))
        elif i%2 == 1:
            south.append(np.fromstring(lines[i], dtype=np.float32, sep=' '))
        else:
            assert(False)

    south = np.array(south)
    north = np.array(north)
    assert(south.shape == (16471, 22))
    assert(north.shape == (16471, 22))

    # note in the SWMF source code, it would appear  Phi <==> Psi  , they are two different names for same thing...
    # likely done at some point to avoid confusion with PHI (the electric potential)
    if True:
        psiN = north[:,4]
        thetaN = north[:,3]
        psiS = south[:,4]
        thetaS = south[:,3]

        psi_ax = 2.*np.arange(181, dtype=np.float32)
        psi2 = np.empty((181,91), dtype=np.float32)
        psi2[:,:] = psi_ax[:,None]

        assert(np.all( psiN == psi2.flatten() ))
        assert(np.all( psiN == psiS ))

        thetaN_ax = np.arange(91, dtype=np.float32)
        thetaN_ax[0] = 5.7296e-3
        thetaN2 = np.empty((181,91), dtype=np.float32)
        thetaN2[:,:] = thetaN_ax[None,:]

        assert(np.all( thetaN == thetaN2.flatten() ))

        thetaS_ax = 90+np.arange(91, dtype=np.float32)
        thetaS_ax[-1] = 179.99
        thetaS2 = np.empty((181,91), dtype=np.float32)
        thetaS2[:,:] = thetaS_ax[None,:]

        assert(np.all( thetaS == thetaS2.flatten() ))

    # note, there are 32942==181*182 points in the data file, but only 2+180*179 are unique, if you dont count the fudging done at the poles
    #   north pole, south pole, and 180x179 PsixTheta grid with Psi <0,2,4,...,356,358> and Theta <1,2,3,...,178,179> 
    # also we want 23 variables, one extra for the integration measure.
    data_arr = np.empty((23, 2+180*179), dtype=np.float32)
    for i in range(22):
        assert(np.all( north[:,i].reshape((181,91))[:, 90] == south[:,i].reshape((181,91))[:, 0] ))
        if i!=4:
            test = np.max(np.abs( north[:,i].reshape((181,91))[0, :] - north[:,i].reshape((181,91))[-1, :] ))
            if test>2e-15:
                print('WARNING ' \
                 +'np.max(np.abs( north[:,i].reshape((181,91))[0, :] - north[:,i].reshape((181,91))[-1, :] )) > 2e-15' \
                 +f'is actually {test}')

            test = np.max(np.abs( south[:,i].reshape((181,91))[0, :] - south[:,i].reshape((181,91))[-1, :] ))
            if test>2e-15:
                print('WARNING ' \
                 +'np.max(np.abs( south[:,i].reshape((181,91))[0, :] - south[:,i].reshape((181,91))[-1, :] )) > 2e-15' \
                 +f'is actually {test}')

        # the following are not the exact same point, due to theta being 0.0057296 and 179.99 instead of 0.0 and 180.0 respectively.
        # But we ignore the fudging, for the purpose of obtaing a good measure, and average over them and asign to a unique point.
        # One for north pole and one for south pole. 
        if i == 3:
            data_arr[i, 0] = 0.
            data_arr[i, 1] = 180.
        elif i==4:
            data_arr[i, 0] = 0. # or np.nan?
            data_arr[i, 1] = 0.
        else:
            data_arr[i, 0] = np.average(north[:,i].reshape((181,91))[:, 0])
            data_arr[i, 1] = np.average(south[:,i].reshape((181,91))[:, -1])

        data_arr[i, 2:2+89*180] = north[:,i].reshape((181,91))[:-1, 1:-1].ravel()
        data_arr[i, 2+89*180:] = south[:,i].reshape((181,91))[:-1, :-1].ravel()

    # from swmf's PostIONO.f90
    Radius = (6378.+100.)/6378.

    deg = np.pi/180.

    if False:
        x_overwrite = Radius*np.cos(data_arr[varidx['Psi'], :]*deg)*np.sin(data_arr[varidx['Theta'], :]*deg)
        y_overwrite = Radius*np.sin(data_arr[varidx['Psi'], :]*deg)*np.sin(data_arr[varidx['Theta'], :]*deg)
        z_overwrite = Radius*np.cos(data_arr[varidx['Theta'], :]*deg)

        data_arr[varidx['X'], :] = x_overwrite
        data_arr[varidx['Y'], :] = y_overwrite
        data_arr[varidx['Z'], :] = z_overwrite

    # integration measure of dA = dTheta*dPhi*(Rad**2)*sin(Theta)
    varidx['measure'] = np.int32(22)
    data_arr[22, :] = (1.*deg)*(2.*deg)*(Radius**2)*np.sin(data_arr[varidx['Theta'], :])

    return data_arr, varidx, units
