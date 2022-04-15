def test_read():


    from os.path import exists
    import numpy as np
    from urllib.request import urlretrieve
    from swmf_file_reader.batsrus_class import get_class_from_native
    from swmf_file_reader.batsrus_class import get_class_from_cdf

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
    batsnative = get_class_from_native(tmpdir + filebase)


    # Get a 513th value of x, y, z, and rho
    var_dict = dict(batsnative.varidx)
    rho = batsnative.data_arr[:, var_dict['rho']][513]
    x = batsnative.data_arr[:,var_dict['x']][513]
    y = batsnative.data_arr[:,var_dict['y']][513]
    z = batsnative.data_arr[:,var_dict['z']][513]
    print(x, y, z, rho)
    # -71.25 -15.75 -7.75 4.22874
    
    # Get interpolated value at a native grid point
    rhoi = batsnative.interpolate(np.array([x, y, z]), 'rho')
    
    assert rho == rhoi

    # Read same data stored in a Kameleon-generated CDF
    filename = '3d__var_2_e20190902-041000-000.out.cdf'
    
    if not exists(tmpdir + filename):
        print("Downloading " + urlbase + filename)
        urlretrieve(urlbase + filename, tmpdir + filename)

    batscdf = get_class_from_cdf(tmpdir + filename)
    
    assert batsnative.data_arr.shape == batscdf.data_arr.shape

