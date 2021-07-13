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

# tested only for DIPTSUR2
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
    if True:
        x_overwrite = Radius*np.cos(data_arr[varidx['Psi'], :]*deg)*np.sin(data_arr[varidx['Theta'], :]*deg)
        y_overwrite = Radius*np.sin(data_arr[varidx['Psi'], :]*deg)*np.sin(data_arr[varidx['Theta'], :]*deg)
        z_overwrite = Radius*np.cos(data_arr[varidx['Theta'], :]*deg)

        data_arr[varidx['X'], :] = x_overwrite
        data_arr[varidx['Y'], :] = y_overwrite
        data_arr[varidx['Z'], :] = z_overwrite

    # integration measure of dA = dTheta*dPhi*(Rad**2)*sin(Theta)
    varidx['measure'] = np.int32(22)
    data_arr[22, :] = (1.*deg)*(2.*deg)*(Radius**2)*np.sin(deg*data_arr[varidx['Theta'], :])

    return data_arr, varidx, units

# tested only for SWPC_SWMF_052811_2
def read_iono_cdf(filename):
    import cdflib

    cdf = cdflib.CDF(filename)

    x = cdf.varget("x")
    y = cdf.varget("y")
    z = cdf.varget("z")
    psi_ax = cdf.varget("phi")
    theta_ax = cdf.varget("theta") # for some reason these seem to be 
    if np.max(theta_ax) < 91:
        theta_ax = 90. - theta_ax # change latitude to colatitude\

    x = x.reshape((181, 181))
    y = y.reshape((181, 181))
    z = z.reshape((181, 181))
    psi_ax = psi_ax.reshape((181,))
    theta_ax = theta_ax.reshape((181,))

    # The following will verify forall (i,j), 
    #       x[i,j] == x_func(psi_ax[i], theta_ax[j])
    #   where x_fuct(psi, theta) is standard cartesian function of polar coordinates
    if True:
        deg = np.pi/180.
        for i in range(181):
            for j in range(181):
                x_calc = np.cos(psi_ax[i]*deg)*np.sin(theta_ax[j]*deg)
                y_calc = np.sin(psi_ax[i]*deg)*np.sin(theta_ax[j]*deg)
                z_calc =                       np.cos(theta_ax[j]*deg)
                #print( np.abs(x[i,j] - x_calc) < 1e-4 )
                #print( np.abs(y[i,j] - y_calc) < 1e-4 )
                #print( np.abs(z[i,j] - z_calc) < 1e-4 )
                #print( (x[i,j] , x_calc) )
                #print( (y[i,j] , y_calc) )
                #print( (z[i,j] , z_calc) )
                assert( np.abs(x[i,j] - x_calc) < 1e-4 )
                assert( np.abs(y[i,j] - y_calc) < 1e-4 )
                assert( np.abs(z[i,j] - z_calc) < 1e-4 )
    # we could therefore in principle build an interpolator
    # for any variable var as a function of polar coordinates using:
    #   var(psi_ax[i], theta_ax[j]) == cdf.varget('var').reshape((181, 181))[i,j]

    #conv = {
    #    'x'                         : 'X'     ,
    #    'y'                         : 'Y'     ,
    #    'z'                         : 'Z'     ,
    #    'theta'                     : 'Theta' ,
    #    'phi'                       : 'Psi'   , # only one that differs by more than just the case
    #    'sigmaH'                    : 'SigmaH',
    #    'sigmaP'                    : 'SigmaP',
    #    'eflux'                     : 'eflux' ,
    #    'eave'                      : 'eave'  ,
    #    'jr'                        : 'Jr'    ,
    #    'ep'                        : 'Ep'    ,
    #    'ex'                        : 'Ex'    ,
    #    'ey'                        : 'Ey'    ,
    #    'ez'                        : 'Ez'    ,
    #    'jx'                        : 'Jx'    ,
    #    'jy'                        : 'Jy'    ,
    #    'jz'                        : 'Jz'    ,
    #    'ux'                        : 'Ux'    ,
    #    'uy'                        : 'Uy'    ,
    #    'uz'                        : 'Uz'    ,
    #    'kameleon_identity_unknown1': None    ,
    #    'kameleon_identity_unknown2': None    ,
    #    }

    varidx = Dict.empty(
        key_type=types.unicode_type,
        value_type=types.int64,
        )
    units = {}

    data_arr = np.empty((23, 2+180*179), dtype=np.float32); data_arr[:,:]=np.nan
    iVar = 0
    for cdfvar in cdf.cdf_info()['zVariables']:
        var = cdf.varattsget(cdfvar)['Original Name']
        units[var] = cdf.varattsget(cdfvar)['units']

        if var=='Psi' or var=='Theta':
            continue
        vararr = cdf.varget(cdfvar).reshape((181, 181)) # note, varget is case insensitive
        varidx[var] = iVar

        data_arr[iVar, 0] = np.average(vararr[i, 180])
        data_arr[iVar, 1] = np.average(vararr[i,   0])

        indx = 2
        for i in range(180):
            for j in range(1,180):
                data_arr[iVar, indx] = vararr[i,j]
                indx += 1
        iVar += 1

    varidx['Psi'] = _Psi = iVar
    varidx['Theta'] = _Theta = iVar+1
    varidx['measure'] = _measure = iVar+2
    assert(_measure == 22)

    data_arr[_Psi, 0] = 0. # or np.nan?
    data_arr[_Psi, 1] = 0.

    data_arr[_Theta, 0] = 0.
    data_arr[_Theta, 1] = 180.

    indx = 2
    for i in range(180):
        for j in range(1,180):
            data_arr[_Psi,   indx] = psi_ax[i]
            data_arr[_Theta, indx] = theta_ax[j]
            indx += 1

    # from swmf's PostIONO.f90
    Radius = (6378.+100.)/6378.

    deg = np.pi/180.
    if True:
        x_overwrite = Radius*np.cos(data_arr[varidx['Psi'], :]*deg)*np.sin(data_arr[varidx['Theta'], :]*deg)
        y_overwrite = Radius*np.sin(data_arr[varidx['Psi'], :]*deg)*np.sin(data_arr[varidx['Theta'], :]*deg)
        z_overwrite = Radius*np.cos(data_arr[varidx['Theta'], :]*deg)
        assert(np.max(np.abs(Radius*data_arr[varidx['X'], :] - x_overwrite))<1e-3)
        assert(np.max(np.abs(Radius*data_arr[varidx['Y'], :] - y_overwrite))<1e-3)
        assert(np.max(np.abs(Radius*data_arr[varidx['Z'], :] - z_overwrite))<1e-3)
        data_arr[varidx['X'], :] = x_overwrite
        data_arr[varidx['Y'], :] = y_overwrite
        data_arr[varidx['Z'], :] = z_overwrite

    # integration measure of dA = dTheta*dPhi*(Rad**2)*sin(Theta)
    data_arr[_measure, :] = (1.*deg)*(2.*deg)*(Radius**2)*np.sin(deg*data_arr[varidx['Theta'], :])

    return data_arr, varidx, units


def testing():
    fname = '/home/gary/temp/i_e20190902-041100-000.tec'
    data_arr, varidx, units, x_overwrite, y_overwrite, z_overwrite = read_iono_tec(fname)
    print(varidx)
    print(units)
    print(data_arr.shape)

    X  = data_arr[varidx['X'] , :]
    Y  = data_arr[varidx['Y'] , :]
    Z  = data_arr[varidx['Z'] , :]
    measure  = data_arr[varidx['measure'] , :]
    Theta  = data_arr[varidx['Theta'] , :]
    Psi  = data_arr[varidx['Psi'] , :]

    R = np.sqrt(X**2 + Y**2 + Z**2)
    print(R)
    # from swmf's PostIONO.f90
    Radius = (6378.+100.)/6378.
    print(Radius)
    print(4*np.pi*Radius**2)
    print(np.sum(measure))

    print(np.max(np.abs(Radius*X-x_overwrite)))
    print(np.max(np.abs(Radius*Y-y_overwrite)))
    print(np.max(np.abs(Radius*Z-z_overwrite)))

    import magnetopost.ionosphere_integrals as mii
    from magnetopost.units_and_constants import phys

    JX = data_arr[varidx['Jx'],:]
    JY = data_arr[varidx['Jy'],:]
    JZ = data_arr[varidx['Jz'],:]

    # note, surface current density (not vol current density)
    KX_p = data_arr[varidx['Ex'],:] * data_arr[varidx['SigmaP'],:]
    KY_p = data_arr[varidx['Ey'],:] * data_arr[varidx['SigmaP'],:]
    KZ_p = data_arr[varidx['Ez'],:] * data_arr[varidx['SigmaP'],:]

    XYZ = data_arr[[varidx['X'],varidx['Y'],varidx['Z']], :].transpose()
    E_iono = data_arr[[varidx['Ex'],varidx['Ey'],varidx['Ez']], :].transpose()
    # note, surface current density (not vol current density)
    unit_b_dipole = mii.get_dipole_field_V(XYZ)
    unit_b_dipole = unit_b_dipole/np.linalg.norm(unit_b_dipole, axis=1)[:,None]
    K_h = np.cross(unit_b_dipole, E_iono) * data_arr[varidx['SigmaH'],:][:,None]
    Measure = data_arr[varidx['measure'],:]

    K_h = K_h   * ( phys['mu0']*phys['Siemens']*phys['mV']/phys['m'] )
    KX_p = KX_p * ( phys['mu0']*phys['Siemens']*phys['mV']/phys['m'] )
    KY_p = KY_p * ( phys['mu0']*phys['Siemens']*phys['mV']/phys['m'] )
    KZ_p = KZ_p * ( phys['mu0']*phys['Siemens']*phys['mV']/phys['m'] )

    JX = JX * ( phys['mu0']*phys['muA']/(phys['m']**1) )
    JY = JY * ( phys['mu0']*phys['muA']/(phys['m']**1) )
    JZ = JZ * ( phys['mu0']*phys['muA']/(phys['m']**1) )

    print('\n\nX\n\n')
    print(K_h[:,0])
    print(KX_p)
    print(-K_h[:,0]+KX_p)
    print(JX)
    print(np.max(np.abs(-K_h[:,0] + KX_p - JX)))


    print('\n\nY\n\n')
    print(K_h[:,1])
    print(KY_p)
    print(-K_h[:,1]+KY_p)
    print(JY)
    print(np.max(np.abs(-K_h[:,1] + KY_p - JY)))

    print('\n\nZ\n\n')
    print(K_h[:,2])
    print(KZ_p)
    print(-K_h[:,2]+KZ_p)
    print(JZ)
    print(np.max(np.abs(-K_h[:,2] + KZ_p - JZ)))

    JREC_X = -K_h[:,0] + KX_p 
    JREC_Y = -K_h[:,1] + KY_p 
    JREC_Z = -K_h[:,2] + KZ_p

    JREC_norm = np.sqrt(JREC_X**2 + JREC_Y**2 + JREC_Z**2)
    J_norm = np.sqrt(JX**2 + JY**2 + JZ**2)

    print(JREC_norm)
    print(J_norm)
    print(np.max(np.abs(J_norm - JREC_norm)))

    import matplotlib.pyplot as plt
    plt.plot(J_norm,'.')
    plt.plot(JREC_norm, '.')
    plt.show()


def testing2():
    fname = '/home/gary/Documents/code_repos/magnetosphere/data/SWPC_SWMF_052811_2/IONO-2D_CDF/SWPC_SWMF_052811_2.swmf.it061214_071000_000.cdf'
    data_arr, varidx, units, x_overwrite, y_overwrite, z_overwrite = read_iono_cdf(fname)
    print(varidx)
    print(units)
    print(data_arr.shape)

    X  = data_arr[varidx['X'] , :]
    Y  = data_arr[varidx['Y'] , :]
    Z  = data_arr[varidx['Z'] , :]
    measure = data_arr[varidx['measure'] , :]
    Theta  = data_arr[varidx['Theta'] , :]
    Psi  = data_arr[varidx['Psi'] , :]

    R = np.sqrt(X**2 + Y**2 + Z**2)
    print(R)
    # from swmf's PostIONO.f90
    Radius = (6378.+100.)/6378.
    print(Radius)
    print('surface area:')
    print(4*np.pi*Radius**2)
    print(np.sum(measure))

    print( np.max(np.abs(Radius*X-x_overwrite)) )
    print( np.max(np.abs(Radius*Y-y_overwrite)) )
    print( np.max(np.abs(Radius*Z-z_overwrite)) )

    print(np.where(np.abs(Radius*X-x_overwrite) > 5e-5))
    print( np.max(np.abs(Radius*X-x_overwrite)[2:]) )
    print( np.max(np.abs(Radius*Y-y_overwrite)[2:]) )
    print( np.max(np.abs(Radius*Z-z_overwrite)[2:]) )
    
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    c = data_arr[varidx['SigmaP'] , :]
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    #sc1 = ax.scatter(X, Y, Z, c=c, s=0.3)
    sc1 = ax.scatter(Radius*X-x_overwrite, Radius*Y-y_overwrite, Radius*Z-z_overwrite, c=c, s=10.)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_xlim((-1.5e-4,1.5e-4))
    ax.set_ylim((-1.5e-4,1.5e-4))
    ax.set_zlim((-1.5e-4,1.5e-4))
    fig.colorbar(sc1, ax=ax, label='SigmaP')
    plt.show()

if __name__ == '__main__':
    filename = '/home/gary/Documents/code_repos/magnetosphere/data/SWPC_SWMF_052811_2/IONO-2D_CDF/SWPC_SWMF_052811_2.swmf.it061214_071000_000.cdf'
    tup=read_iono_cdf(filename)
    print(tup[0])
    print(tup[1])
    print(tup[2])
    #testing2()
