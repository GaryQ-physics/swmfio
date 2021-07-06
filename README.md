# swmf_file_reader

reads the .out, .tree, and .info files for the data from the
(BATSRUS) magnetosphere component of an SWMF run.
For example datafiles, go to
`http://mag.gmu.edu/git-data/GaryQ-Physics/demodata/`.

Full functionality for python3, limited functionality python2 due to no
python2 compatibility with numba.

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

# using this package

## Batsrus files

### overview

Bats-r-us outputs consists of files: 3d_*.out , 3d_*.tree , 3d_*.info

3d_*.out:
> consists of a sequence of 1d arrays of the same length and corresponding variable names,
> together with a few extra parameter.
> Those variable names include x, y, and z coordinates.
> Each entry of the arrays corresponds to a unique gridpoint
> (and has a unique (x,y,z) value).

3d_*.tree:
> Consists of a 2d arrays describing the block tree, and a few extra parameters.
> The 2d array is refered to internally as "iTree_IA" (following the BATRUS source code name).
> The 2d array describes a tree, and the nodes of the tree are indexed by "iNode".
> Not all the nodes are used to get the gridpoints, only the leaves of the tree are.
> Restricting to just the leaves  yield "iNodes" that are not a simple range (1,2,3,...),
> so they can alternatively be indexed by a different "iBlock".
> Each leaf corresponds to a nI-by-nJ-by-nK regular grid of gridpoints,
> where "nI","nJ",and "nK" are fixed and determined from the other files.

3d_*.info:
> text file containing nI,nJ,nK, and other meta data.



There are two interfaces for working with theses files:
batsrus_imperative and batsrus_class

### class

```
from swmf_file_reader.batsrus_class import return_class
batsclass = return_class('/tmp/3d__var_2_e20190902-041000-000')
```

This returns a numba jit class,
with the simulation data in arrays as class atributes,
and interpolation and differentiation as class methods

### imperative

```
from swmf_file_reader.batsrus_imperative import read_files
batstup = read_files('/tmp/3d__var_2_e20190902-041000-000')
```

This returns a python Named Tuple,
containing the simulation data in arrays ect.
To do interpolation and differentiation, you must pass the namedtuple
to the associated function is batsrus_imperative.
NOTE: this can be used with python2, although it will be slower since
numba will be bypassed.

### Basic Example

```
filetag = '/tmp/3d__var_2_e20190902-041000-000'

########## class
from swmf_file_reader.batsrus_class import return_class
batsclass = return_class(filetag)

print( batsclass.data_arr.shape )
print( batsclass.interpolate(np.array([1.,1.,1.), 'rho') )
print( batsclass.get_native_partial_derivatives(123456, 'rho') )

########## imperative
from swmf_file_reader.batsrus_imperative import read_files, get_native_partial_derivatives, interpolate
batstup = read_files(filetag)

print( batstup.data_arr.shape )
print( interpolate(batstup, np.array([1.,1.,1.), 'rho') )
print( get_native_partial_derivatives(batstup, 123456, 'rho') )
```

### export vtk files

```
from swmf_file_reader.swmf2vtk import write_BATSRUS_unstructured_grid_vtk
filetag = '/tmp/3d__var_2_e20190902-041000-000'
write_BATSRUS_unstructured_grid_vtk(filetag, use_ascii=False)
```
This will create '/tmp/3d__var_2_e20190902-041000-000.vtk'
Note, currently uses the batsrus_class, so only for python3.
