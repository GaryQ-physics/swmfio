import numpy as np
from swmf_file_reader.read_ie_files import read_iono_cdf, read_iono_tec
from swmf_file_reader.batsrus_class import get_class_from_native


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


def testing_cdf():
    import cdflib
    cdf = cdflib.CDF('/home/gary/Documents/code_repos/magnetosphere/data/SWPC_SWMF_052811_2/GM_CDF/3d__var_1_t00001001_n0002710.out.cdf')

    assert(1968*8**3 == 1007616)
    assert(1968 == cdf.varget('block_at_amr_level')[0,0])
    assert(1968 == np.where(cdf.varget('block_amr_levels') == 0)[1][0])
    assert(1968 == cdf.globalattsget()['number_of_blocks'])

    for iNodeP in range(cdf.varget('block_amr_levels').size):
        if iNodeP < 1968:
            assert(cdf.varget('block_child_count')[0,iNodeP] == 0)
        else:
            assert(cdf.varget('block_child_count')[0,iNodeP] == 8)

    print(cdf.varget('block_amr_levels')[0,1968])
    
    print(cdf.varget('block_child_id_1')[0,1968])
    print(cdf.varget('block_child_id_2')[0,1968])
    print(cdf.varget('block_child_id_3')[0,1968])
    print(cdf.varget('block_child_id_4')[0,1968])
    print(cdf.varget('block_child_id_5')[0,1968])
    print(cdf.varget('block_child_id_6')[0,1968])
    print(cdf.varget('block_child_id_7')[0,1968])
    print(cdf.varget('block_child_id_8')[0,1968])
    
    print(cdf.varget('block_amr_levels')[0,1969])
    print(cdf.varget('block_amr_levels')[0,1970])
    print(cdf.varget('block_amr_levels')[0,1971])
    print(cdf.varget('block_amr_levels')[0,1972])
    print(cdf.varget('block_amr_levels')[0,1973])
    print(cdf.varget('block_amr_levels')[0,1974])
    print(cdf.varget('block_amr_levels')[0,1975])
    print(cdf.varget('block_amr_levels')[0,1976])

    print(cdf.varget('block_x_min')[0,1968])
    print(cdf.varget('block_y_min')[0,1968])
    print(cdf.varget('block_z_min')[0,1968])
    print(cdf.varget('block_x_max')[0,1968])
    print(cdf.varget('block_y_max')[0,1968])
    print(cdf.varget('block_z_max')[0,1968])


    print(cdf.varget('block_x_min')[0,1969])
    print(cdf.varget('block_y_min')[0,1969])
    print(cdf.varget('block_z_min')[0,1969])
    print(cdf.varget('block_x_max')[0,1969])
    print(cdf.varget('block_y_max')[0,1969])
    print(cdf.varget('block_z_max')[0,1969])


def test():
    #cls = get_class_from_cdf('/home/gary/Documents/code_repos/magnetosphere/data/SWPC_SWMF_052811_2/GM_CDF/3d__var_1_t00001001_n0002710.out.cdf')
    cls = get_class_from_native('/tmp/3d__var_2_e20190902-041000-000')
    print(cls.varidx)
    print(cls.data_arr)

    onpoint = np.array([-146.,  -14.,  -14.])
    offpoint = np.array([-143.,  -15.,  -15.])
    print(cls.interpolate(onpoint, 'rho'))
    print(cls.interpolate(onpoint, 'x'))
    print(cls.interpolate(onpoint, 'y'))
    print(cls.interpolate(onpoint, 'z'))
    print(cls.interpolate(offpoint, 'rho'))
    print(cls.interpolate(offpoint, 'x'))
    print(cls.interpolate(offpoint, 'y'))
    print(cls.interpolate(offpoint, 'z'))
    print(cls.get_native_partial_derivatives(123456, 'rho'))
    print(cls.get_native_partial_derivatives(125356, 'rho'))
    print(cls.get_native_partial_derivatives(143456, 'rho'))
    print(cls.get_native_partial_derivatives(143456, 'x'))
    print(cls.get_native_partial_derivatives(143456, 'y'))
    print(cls.get_native_partial_derivatives(143456, 'z'))
