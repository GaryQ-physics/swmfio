# Overview

Reads the `.out`, `.tree`, and `.info` files for the data from the (BATSRUS) magnetosphere component of an SWMF run and provides interpolator method.

Can also read CCMC `.cdf`, which contain the most of the information contained in the `.out`, `.tree`, and `.info` files.

Also provides a function to output data on native grid (as unstructured grid with native connectivity) to a VTK file.

For example data files, see `http://mag.gmu.edu/git-data/GaryQ-Physics/demodata/`.

Requires Python 3.

# Install

for user:
```
pip install 'git+https://github.com/GaryQ-physics/swmf_file_reader.git' --upgrade
```

for developer:
```
git clone https://github.com/GaryQ-physics/swmf_file_reader.git
cd swmf_file_reader
pip install --editable .
```

# Using

## BATSRUS files

### Overview

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

### Usage

```
from swmf_file_reader.batsrus_class import return_class
batsclass = return_class('/tmp/3d__var_2_e20190902-041000-000')
```

This returns a numba jit class, with the simulation data in arrays as class atributes, and interpolation and differentiation as class methods

### Example

```
filetag = '/tmp/3d__var_2_e20190902-041000-000'

from swmf_file_reader.batsrus_class import return_class
batsclass = return_class(filetag)

print( batsclass.data_arr.shape )
print( batsclass.interpolate(np.array([1.,1.,1.), 'rho') )
print( batsclass.get_native_partial_derivatives(123456, 'rho') )

## CDF files


## Export VTK file

```
from swmf_file_reader.swmf2vtk import write_BATSRUS_unstructured_grid_vtk
filetag = '/tmp/3d__var_2_e20190902-041000-000'
write_BATSRUS_unstructured_grid_vtk(filetag, use_ascii=False)
```

This will create `/tmp/3d__var_2_e20190902-041000-000.vtk`
