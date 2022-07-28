import swmfio

import logging
swmfio.logger.setLevel(logging.INFO)

demo_num = 2

if demo_num == 0:
    url = 'http://mag.gmu.edu/git-data/swmfio/3d__var_1_t00000000_n0002500.out.cdf'
    # Downloads .cdf file if local copy not found
    file = swmfio.dlfile(url, progress=True)
    batsclass = swmfio.read_batsrus(file)
    vtkfile = swmfio.write_vtk(batsclass)

if demo_num == 1:
    url = 'http://mag.gmu.edu/git-data/swmfio/3d__var_1_t00000000_n0002500.out.cdf'
    # Downloads .cdf file if local copy not found
    file = swmfio.dlfile(url, progress=True)
    vtkfile = swmfio.write_vtk(file)
    print(vtkfile) # /tmp/mag.gmu.edu/git-data/swmfio/3d__var_1_t00000000_n0002500.out.cdf.vtk

if demo_num == 2:
    # .out, .tree, and .info file are remote.
    url = 'http://mag.gmu.edu/git-data/swmfio/3d__var_2_e20190902-041000-000'
    # Downloads .out, .tree, and .info files if local copy not found
    filebase = swmfio.dlfile(url, progress=True)
    vtkfile = swmfio.write_vtk(filebase)
    print(vtkfile) # /tmp/mag.gmu.edu/git-data/swmfio/d__var_2_e20190902-041000-000.vtk

if demo_num == 3:
    # .out, .tree, and .info file are local
    filebase = '/tmp/3d__var_2_e20190902-041000-000'
    vtkfile = swmfio.write_vtk(filebase)
    print(vtkfile) # /tmp/3d__var_2_e20190902-041000-000.vtk

if demo_num == 4:
    # .out, .tree, and .info file are local
    filebase = '/tmp/3d__var_2_e20190902-041000-000'
    vtkfile = swmfio.write_vtk(filebase, variables=['b'])
    print(vtkfile) # /tmp/3d__var_2_e20190902-041000-000_vars=b.vtk

if demo_num == 5:
    # For debugging, can output data for selected blocks and output ASCII VTK.
    swmfio.logger.setLevel(logging.DEBUG)
    filebase = '/tmp/3d__var_2_e20190902-041000-000'
    vtkfile = swmfio.write_vtk(filebase, use_ascii=True, blocks=[0, 1])
    print(vtkfile) # /tmp/3d__var_2_e20190902-041000-000_blocks=0,1.vtk
