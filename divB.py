import swmfio

import logging
swmfio.logger.setLevel(logging.INFO)

# swmfio binary reader fails
#url = 'http://mag.gmu.edu/git-data/dwelling/divB_simple1/GM/3d__mhd_4_e20100320-013000-000'

# works
#url = 'http://mag.gmu.edu/git-data/swmfio/3d__var_2_e20190902-041000-000'

# infinite loop in find_tree_node()
url = 'http://mag.gmu.edu/git-data/bcurtiswx/Differences/data/Brian_Curtis_042213_2/GM_CDF/3d__var_1_e20000101-002500-000.out.cdf'
file = swmfio.dlfile(url, progress=True)
batsclass = swmfio.read_batsrus(file)
vtkfile = swmfio.write_vtk(batsclass)
