import swmfio

url = 'http://mag.gmu.edu/git-data/swmfio/i_e20190902-041100-000.tec'

file = swmfio.dlfile(url, progress=True)

data_arr, varidx, units = swmfio.read_rim(file)

print("swmfio.read_rim of tec file")
print(data_arr)
print(varidx)
print(units)

url = 'http://mag.gmu.edu/git-data/swmfio/it100320_000100_000.idl'

file = swmfio.dlfile(url, progress=True)

from spacepy.pybats import rim
pbo = rim.Iono(file)
pbo.readascii()

names = pbo.keys()
units = pbo.listunits()

print("swmfio.read_rim of idl file")
print(names)
print(units)
#data_arr = np.empty((23, 2+180*179), dtype=np.float32)
x = pbo.get('s_jouleheat')
print(x.shape)