import swmfio

import logging
swmfio.logger.setLevel(logging.INFO)

if True:
    # If .out, .tree, and .info file are remote.
    url = 'http://mag.gmu.edu/git-data/swmfio/3d__var_2_e20190902-041000-000'
    # Downloads .out, .tree, and .info files if not found there already.
    vtkfile = swmfio.write_vtk(url)
    print(vtkfile) # /tmp/mag.gmu.edu/git-data/swmfio/d__var_2_e20190902-041000-000.vtk

if False:
    # If .out, .tree, and .info file are local
    file = '/tmp/3d__var_2_e20190902-041000-000'
    vtkfile = swmfio.write_vtk(file)
    print(vtkfile) # /tmp/3d__var_2_e20190902-041000-000.vtk

if False:
    # If .out, .tree, and .info file are local
    file = '/tmp/3d__var_2_e20190902-041000-000'
    vtkfile = swmfio.write_vtk(file, variables=['b'])
    print(vtkfile) # /tmp/3d__var_2_e20190902-041000-000_vars=b.vtk


if False:
    # For debugging, can output data for selected blocks
    # and output ASCII file.
    swmfio.logger.setLevel(logging.DEBUG)
    file = '/tmp/3d__var_2_e20190902-041000-000'
    vtkfile = swmfio.write_vtk(file, use_ascii=True, blocks=[0, 1])
    print(vtkfile) # /tmp/3d__var_2_e20190902-041000-000_blocks=0,1.vtk
