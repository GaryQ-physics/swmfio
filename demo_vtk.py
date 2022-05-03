from os.path import exists
from urllib.request import urlretrieve

import swmf_file_reader as swmf

urlbase = 'http://mag.gmu.edu/git-data/swmf_file_reader/demodata/'
tmpdir = '/tmp/'
filebase = '3d__var_2_e20190902-041000-000'

for ext in ['.tree', '.info', '.out']:
    filename = tmpdir + filebase + ext
    if not exists(filename):
        print("Downloading " + urlbase + filename)
        print("to")
        print(tmpdir + filename)
        urlretrieve(urlbase + filename, tmpdir + filename)

swmf.write_vtk(tmpdir + filebase, debug=True)
