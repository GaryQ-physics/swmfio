import numpy as np
from urllib.request import urlretrieve
from swmf_file_reader.batsrus_class import get_class_from_cdf

urlbase = 'http://mag.gmu.edu/git-data/swmf_file_reader/demodata/'
tmpdir = '/tmp/'
filename = '3d__var_1_t00000000_n0002500.out.cdf'

print("Downloading " + urlbase + filename)
urlretrieve(urlbase + filename, tmpdir + filename)

batsclass = get_class_from_cdf(tmpdir + filename)

print( batsclass.data_arr.shape )
# (1007616, 20)

print( batsclass.interpolate(np.array([1.,1.,1.]), 'rho') )
# 28.0

print( batsclass.get_native_partial_derivatives(123456, 'rho') )
# [ 1.8113604  -0.00352001  2.1863403 ]