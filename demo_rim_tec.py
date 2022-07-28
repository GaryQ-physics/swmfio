import swmfio as swmfio
from os.path import exists
from urllib.request import urlretrieve

url = 'http://mag.gmu.edu/git-data/swmfio/i_e20190902-041100-000.tec'

file = swmfio.dlfile(url, progress=True)

data, varidx, units = swmfio.read_rim(file)

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