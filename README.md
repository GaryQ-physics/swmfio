# Overview

Reads some magnetosphere and ionosphere data files from an SWMF run.
For the (BATSRUS) magnetosphere module, reads `.out`, `.tree`, and `.info` files.
For the (RIM) ionosphere module, reads `.tec` files.

Can also read CCMC `.cdf` files, which contain the most the information contained in the native magnetosphere files. CDF files are read using [cdflib](`https://pypi.org/project/cdflib/`), which is not installed by default.

For BATSRUS: returns a Numba class, with the simulation data in arrays as class atributes,
and interpolation and differentiation as class methods.
For RIM: returns a tuple `(data_arr, varidx, units)`,
where `data_arr` is a numpy array with the data,
`varidx` is a Numba typed Dictionary which maps variable name strings to their corresponding index in `data_arr`,
and `units` is a Dictionary of the corresponding units.

This code also provides a function to output the magnetosphere (BATSRUS) data on native grid to a VTK file, as either an unstructured voxel or hexahedra grid.

This package has a similar functionality to the Julia package [Batsrus.jl](https://github.com/henry2004y/Batsrus.jl)
with the exceptions that in `swmf_file_reader`
* the VTK grid for  places the cell-centered data at the center of cells instead of at cell vertices.
* there is an interpolator interface, and
* CCMC `.cdf` files can be read.

For example data files, see [http://mag.gmu.edu/git-data/swmf_file_reader/demodata/](http://mag.gmu.edu/git-data/swmf_file_reader/demodata/).

This code is used in [https://github.com/GaryQ-physics/magnetopost](https://github.com/GaryQ-physics/magnetopost) to post process magnetosphere simulation data.

# Install

Requires Python 3.

## User:

```
pip install 'git+https://github.com/GaryQ-physics/swmf_file_reader.git' --upgrade
# pip install cdfilb 
```

## Developer:
```
git clone https://github.com/GaryQ-physics/swmf_file_reader.git
cd swmf_file_reader
pip install --editable .
# pip install cdfilb 
```

# Background

BATSRUS outputs consists of files: `3d_*.out`, `3d_*.tree`, `3d_*.info`

`3d_*.out:`
> consists of a sequence of 1d arrays of the same length and corresponding variable names,
> together with a few extra parameter.
> Those variable names include x, y, and z coordinates.
> Each entry of the arrays corresponds to a unique gridpoint
> (and has a unique (x,y,z) value).

`3d_*.tree:`
> Consists of a 2d arrays describing the block tree, and a few extra parameters.
> The 2d array is refered to internally as "iTree_IA" (following the BATRUS source code name).
> The 2d array describes a tree, and the nodes of the tree are indexed by "iNode".
> Not all the nodes are used to get the gridpoints, only the leaves of the tree are.
> Restricting to just the leaves  yield "iNodes" that are not a simple range (1,2,3,...),
> so they can alternatively be indexed by a different "iBlock".
> Each leaf corresponds to a nI-by-nJ-by-nK regular grid of gridpoints,
> where "nI","nJ",and "nK" are fixed and determined from the other files.

`3d_*.info:`
> text file containing nI,nJ,nK, and other meta data.

RIM output consists of file: `i_*.tec`
> A text file, consisting of header with variable information,
> and then the data in ASCII format.
> The data consists of the value of the variables given on a 2d regular grid in latitude longitude.
> The radius is fixed.

# Usage

## BATSRUS datafiles example

Download the demo files by running following in your shell.

```
wget -P /tmp http://mag.gmu.edu/git-data/swmf_file_reader/demodata/3d__var_2_e20190902-041000-000.info
wget -P /tmp http://mag.gmu.edu/git-data/swmf_file_reader/demodata/3d__var_2_e20190902-041000-000.tree
wget -P /tmp http://mag.gmu.edu/git-data/swmf_file_reader/demodata/3d__var_2_e20190902-041000-000.out
```

then run the following in Python.

```python
filetag = '/tmp/3d__var_2_e20190902-041000-000'

from swmf_file_reader.batsrus_class import get_class_from_native
batsclass = get_class_from_native(filetag)

print( batsclass.data_arr.shape )
print( batsclass.interpolate(np.array([1.,1.,1.]), 'rho') )
print( batsclass.get_native_partial_derivatives(123456, 'rho') )
```

## Magnetosphere CDF file example

Download the demo file

```
wget -P /tmp http://mag.gmu.edu/git-data/swmf_file_reader/demodata/3d__var_1_t00000000_n0002500.out.cdf
```

```python
filename = '/tmp/3d__var_1_t00000000_n0002500.out.cdf'

from swmf_file_reader.batsrus_class import get_class_from_cdf
batsclass = get_class_from_cdf(filename)

print( batsclass.data_arr.shape )
print( batsclass.interpolate(np.array([1.,1.,1.]), 'rho') )
print( batsclass.get_native_partial_derivatives(123456, 'rho') )
```

## RIM datafile example

Download the demo file

```
wget  -P /tmp http://mag.gmu.edu/git-data/swmf_file_reader/demodata/i_e20190902-041100-000.tec
```

```python
filename = '/tmp/i_e20190902-041100-000.tec'

from swmf_file_reader.read_ie_files import read_iono_tec
data_arr, varidx, units = read_iono_tec(filename)

print(data_arr.shape)
print(varidx['SigmaP']) # Pedersen conductance

print(data_arr[varidx['Theta'],:]) # the colatitudes
print(data_arr[varidx['Psi'],:]) # the longitudes
```

## Ionosphere CDF file example

Download the demo file

```
wget  -P /tmp http://mag.gmu.edu/git-data/swmf_file_reader/demodata/SWPC_SWMF_052811_2.swmf.it061214_071000_000.cdf
```

```python
filename = '/tmp/SWPC_SWMF_052811_2.swmf.it061214_071000_000.cdf'

from swmf_file_reader.read_ie_files import read_iono_cdf
data_arr, varidx, units = read_iono_cdf(filename)

print(data_arr.shape)
print(varidx['SigmaP']) # Pedersen conductance

print(data_arr[varidx['Theta'],:]) # the colatitudes
print(data_arr[varidx['Psi'],:]) # the longitudes
```

## Export VTK file

Here we already have files
`/tmp/3d__var_2_e20190902-041000-000.out`,
`/tmp/3d__var_2_e20190902-041000-000.tree`,
`/tmp/3d__var_2_e20190902-041000-000.info`.

To do this, you will need magnetovis installed, see `https://github.com/rweigel/magnetovis/`

Run in python:
```
from swmf_file_reader.swmf2vtk import write_BATSRUS_unstructured_grid_vtk
filetag = '/tmp/3d__var_2_e20190902-041000-000'
write_BATSRUS_unstructured_grid_vtk(filetag, use_ascii=False)
```

This will create `/tmp/3d__var_2_e20190902-041000-000.vtk`
