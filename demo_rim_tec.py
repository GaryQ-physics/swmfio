import swmfio as swmfio
from os.path import exists
from urllib.request import urlretrieve

urlbase = 'http://mag.gmu.edu/git-data/swmfio/'
tmpdir = '/tmp/'
filename = 'i_e20190902-041100-000.tec'

if not exists(tmpdir + filename):
    print("Downloading " + urlbase + filename)
    print("to")
    print(tmpdir + filename)
    urlretrieve(urlbase + filename, tmpdir + filename)

data, varidx, units = swmfio.read_rim(tmpdir + filename)

import pprint; pp = pprint.PrettyPrinter(indent=4)
pp.pprint(dict(varidx))
"""
{   'Ave-E': 8,
    'E-Flux': 7,
    'Ex': 11,
    'Ey': 12,
    'Ez': 13,
    'IonNumFlux': 21,
    'JR': 9,
    'JouleHeat': 20,
    'Jx': 14,
    'Jy': 15,
    'Jz': 16,
    'PHI': 10,
    'Psi': 4,
    'SigmaH': 5,
    'SigmaP': 6,
    'Theta': 3,
    'Ux': 17,
    'Uy': 18,
    'Uz': 19,
    'X': 0,
    'Y': 1,
    'Z': 2,
    'measure': 22}
"""

print(varidx['SigmaP']) # Pedersen conductance
# 6
print(varidx['SigmaH']) # Hall conductance
# 5

print(data[varidx['Theta'],:]) # colatitudes
print(data[varidx['Psi'],:])   # longitudes

# [  0. 180.   1. ... 177. 178. 179.]
# [  0.   0.   0. ... 358. 358. 358.]