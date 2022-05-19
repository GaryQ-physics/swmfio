import swmfio as swmfio
from os.path import exists
from urllib.request import urlretrieve

urlbase = 'http://mag.gmu.edu/git-data/swmfio/'
tmpdir = '/tmp/'
filename = 'SWPC_SWMF_052811_2.swmf.it061214_071000_000.cdf'

if not exists(tmpdir + filename):
    print("Downloading " + urlbase + filename)
    print("to")
    print(tmpdir + filename)
    urlretrieve(urlbase + filename, tmpdir + filename)

data, varidx, units = swmfio.read_rim(tmpdir + filename)

import pprint; pp = pprint.PrettyPrinter(indent=4)
pp.pprint(dict(varidx))
"""
{   'Ave-E': 6,
    'E-Flux': 5,
    'Ex': 9,
    'Ey': 10,
    'Ez': 11,
    'IonNumFlux': 19,
    'JR': 7,
    'JouleHeat': 18,
    'Jx': 12,
    'Jy': 13,
    'Jz': 14,
    'PHI': 8,
    'Psi': 20,
    'SigmaH': 3,
    'SigmaP': 4,
    'Theta': 21,
    'Ux': 15,
    'Uy': 16,
    'Uz': 17,
    'X': 0,
    'Y': 1,
    'Z': 2,
    'measure': 22}
"""

print(varidx['SigmaP']) # Pedersen conductance
# 4
print(varidx['SigmaH']) # Hall conductance
# 3

print(data[varidx['Theta'],:]) # colatitudes
#[  0. 180. 179. ...   3.   2.   1.]
print(data[varidx['Psi'],:])   # longitudes
#[  0.   0.   0. ... 358. 358. 358.]