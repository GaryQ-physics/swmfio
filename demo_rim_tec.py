import swmf_file_reader as swmf
from os.path import exists
from urllib.request import urlretrieve

urlbase = 'http://mag.gmu.edu/git-data/swmf_file_reader/demodata/'
tmpdir = '/tmp/'
filename = 'i_e20190902-041100-000.tec'

if not exists(tmpdir + filename):
    print("Downloading " + urlbase + filename)
    print("to")
    print(tmpdir + filename)
    urlretrieve(urlbase + filename, tmpdir + filename)

data_arr, varidx, units = swmf.read_rim(tmpdir + filename)

print(varidx)
# {X: 0, Y: 1, Z: 2, 
#  Theta: 3, Psi: 4, SigmaH: 5, SigmaP: 6, 
#  E-Flux: 7, Ave-E: 8, JR: 9, PHI: 10,
#  Ex: 11, Ey: 12, Ez: 13,
#  Jx: 14, Jy: 15, Jz: 16,
#  Ux: 17, Uy: 18, Uz: 19,
#  JouleHeat: 20, IonNumFlux: 21,
#  measure: 22}

print(varidx['SigmaP']) # Pedersen conductance
print(varidx['SigmaH']) # Hall conductance

print(data_arr[varidx['Theta'],:]) # colatitudes
print(data_arr[varidx['Psi'],:])   # longitudes