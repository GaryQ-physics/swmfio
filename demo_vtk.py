import swmfio

import logging
swmfio.logger.setLevel(logging.INFO)

demo_nums = [1, 2, 3, 4, 5, 6, 7, 8]
demo_nums = [8]

for demo_num in demo_nums:
    if demo_num == 1:
        url = 'http://mag.gmu.edu/git-data/swmfio/3d__var_1_t00000000_n0002500.out.cdf'
        # Downloads .cdf file if local copy not found
        file = swmfio.dlfile(url, progress=True)
        batsclass = swmfio.read_batsrus(file)
        vtkfile = swmfio.write_vtk(batsclass)

    if demo_num == 2:
        url = 'http://mag.gmu.edu/git-data/swmfio/3d__var_1_t00000000_n0002500.out.cdf'
        # Downloads .cdf file if local copy not found
        file = swmfio.dlfile(url, progress=True)
        vtkfile = swmfio.write_vtk(file)
        print(vtkfile) # /tmp/mag.gmu.edu/git-data/swmfio/3d__var_1_t00000000_n0002500.out.cdf.vtk

    if demo_num == 3:
        # .out, .tree, and .info file are remote.
        url = 'http://mag.gmu.edu/git-data/swmfio/3d__var_2_e20190902-041000-000'
        # Downloads .out, .tree, and .info files if local copy not found
        filebase = swmfio.dlfile(url, progress=True)
        vtkfile = swmfio.write_vtk(filebase)
        print(vtkfile) # /tmp/mag.gmu.edu/git-data/swmfio/3d__var_2_e20190902-041000-000.vtk

    if demo_num == 4:
        # .out, .tree, and .info file are local (must run demo 2 first)
        filebase = '/tmp/mag.gmu.edu/git-data/swmfio/3d__var_2_e20190902-041000-000'
        vtkfile = swmfio.write_vtk(filebase)
        print(vtkfile) # /tmp/mag.gmu.edu/git-data/swmfio/3d__var_2_e20190902-041000-000_blocks=0,1.vtk

    if demo_num == 5:
        # .out, .tree, and .info file are local
        filebase = '/tmp/mag.gmu.edu/git-data/swmfio/3d__var_2_e20190902-041000-000'
        vtkfile = swmfio.write_vtk(filebase, variables=['b'])
        print(vtkfile) # /tmp/mag.gmu.edu/git-data/swmfio/3d__var_2_e20190902-041000-000_vars=b.vtk

    if demo_num == 6:
        # For debugging, can output data for selected blocks and output ASCII VTK.
        swmfio.logger.setLevel(logging.DEBUG)
        filebase = '/tmp/mag.gmu.edu/git-data/swmfio/3d__var_2_e20190902-041000-000'
        vtkfile = swmfio.write_vtk(filebase, use_ascii=True, blocks=[0, 1])
        print(vtkfile) # /tmp/mag.gmu.edu/git-data/swmfio/3d__var_2_e20190902-041000-000_blocks=0,1.vtk
        swmfio.logger.setLevel(logging.INFO)

    if demo_num == 7:
        # This file as a tree structure with 3 top-level nodes instead of the usual of 1.
        url = 'http://mag.gmu.edu/git-data/bcurtiswx/Differences/data/Brian_Curtis_042213_2/GM_CDF/3d__var_1_e20000101-002500-000.out.cdf'
        file = swmfio.dlfile(url, progress=True)
        batsclass = swmfio.read_batsrus(file)
        vtkfile = swmfio.write_vtk(batsclass)

    if demo_num == 8:
        url = 'http://mag.gmu.edu/git-data/dwelling/divB_simple1/GM/3d__mhd_4_e20100320-000000-000'
        file = swmfio.dlfile(url, progress=True)
        batsclass = swmfio.read_batsrus(file)
        vtkfile = swmfio.write_vtk(batsclass)
