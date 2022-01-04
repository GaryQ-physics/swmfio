from os.path import exists
from urllib.request import urlretrieve

from swmf_file_reader.swmf2vtk import write_BATSRUS_unstructured_grid_vtk

urlbase = 'http://mag.gmu.edu/git-data/swmf_file_reader/demodata/'
tmpdir = '/tmp/'
filebase = '3d__var_2_e20190902-041000-000'

for ext in ['.tree', '.info', '.out']:
    filename = tmpdir + filebase + ext
    if not exists(filename):
        print("Downloading " + urlbase + filebase + ext)
        urlretrieve(urlbase + filebase + ext, filename)

filebase = tmpdir + filebase
print("Reading {}.*".format(filebase))
write_BATSRUS_unstructured_grid_vtk(filebase)