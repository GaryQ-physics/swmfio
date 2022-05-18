from os.path import exists
from urllib.request import urlretrieve

import swmfio as swmf

urlbase = 'http://mag.gmu.edu/git-data/swmfio/'
tmpdir = '/tmp/'
filebase = '3d__var_2_e20190902-041000-000'

for ext in ['.tree', '.info', '.out']:
    filename = filebase + ext
    if not exists(tmpdir + filename):
        print("Downloading " + urlbase + filename)
        print("to")
        print(tmpdir + filename)
        urlretrieve(urlbase + filename, tmpdir + filename)

# Optional
import logging
swmf.logger.setLevel(logging.INFO)

swmf.write_vtk(tmpdir + filebase, logger=swmf.logger)

if False:
    # For debugging, can output data for selected blocks and output ASCII file.
    import logging
    swmf.logger.setLevel(logging.DEBUG)
    swmf.write_vtk(tmpdir + filebase, logger=swmf.logger, use_ascii=True, blocks=[0, 1])

