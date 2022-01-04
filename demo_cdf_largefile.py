import numpy as np
from urllib.request import urlretrieve
from swmf_file_reader.batsrus_class import get_class_from_cdf

urlbase = 'http://mag.gmu.edu/git-data/swmf_file_reader/demodata/'
tmpdir = '/tmp/'
filename = '3d__var_2_e20190902-041000-000.out.cdf'
    
print("Downloading " + urlbase + filename)
urlretrieve(urlbase + filename, tmpdir + filename)

batsclass = get_class_from_cdf(tmpdir + filename)

print( batsclass.data_arr.shape )
# (5896192, 19)

print( batsclass.interpolate(np.array([1.,1.,1.]), 'rho') )
# 974.5159912109375

print( batsclass.get_native_partial_derivatives(123456, 'rho') )
# [0.25239944 0.41480255 0.7658005 ]