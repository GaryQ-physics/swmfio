# swmf_file_reader

reads the .out, .tree, and .info files for the data from the
(BATSRUS) magnetosphere component of an SWMF run.

In principle works for python 2 and 3

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

# Tested

works:
Python 3.7.9 \[GCC 7.3.0\]
Python 2.7.18 \[GCC 7.3.0\]

# Overview:

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

This package has 4 functions that one would typically use externally,
all located in the 'read_swmf_files.py' file.
The rest are is mostly meant for internal use.
The 4 functions are:
 - read_all(filetag)
    - takes as input the string for 3d_*  without the extentions (including full path if needed).
    - returns a dictionary (refered to internally as "cache")
       which has all the data needed from the 3d_*.tree/out/info files.
       In particular is has the following important key:values
        - "DataArray" : a 5d array.
          the first index corresponds to a variable name from
          the variables in the 3d_*.out file
          (the correspondance is in 'named_var_indexes.py' file).
          After that, the next 4 indexes iterate over the native grid block by block.
          e.g: `DataArray[_rho, iBlockP, i,j,k]` will give you rho at a certain gridpoint.
          That grid point is located in the (i,j,k)th location of the nI-by-nJ-by-nK regular grid
          of the block indexed by iBlock.
          Note: a sublty is the index "iBlockP" starts at 0 (P stands for python)
          and is thus shifted by 1 from "iBlock" (which originally was fortran indexed).
        - "block2node" : a 1d array for converting iBlock indexes to the corresponding iNode indexes.
          The array is indexed by iBlockP, and the entries are the corresponding iNodeP
        - "node2block" : a 1d array for converting iNode indexex to the corresponding iBlock indexes,
          iff iNode is a valid leaf.
          The array is indexed by iNodeP, and the entries are the corresponding iBlockP
          if iNode is a valid leaf, and -1 otherwise.

 - find_index(filetag, point, cache=None, debug=False)
    * Inputed: a (3,) array with the x,y,z of a point in space, a filetag 3d_* and optionally its cache
    * returns the `(iBlockP,i,i,k)` index of point, if point is a gridpoint.
      If point is *not* a gridpoint, then return None.
    * Note: If cache is None, then read_all(filetag) is called,
      but if you already have cache saved you can optionally pass it in
      to save time by not re-loading the datafiles.

 - interpolate(filetag, point, var='p', cache=None, debug=False)
    * Inputed: a (3,) array with the x,y,z of a point in space, a filetag 3d_* and optionally its cache
    * returns the interpolated value of the variable "var" at point.
    * Note: If cache is None, then read_all(filetag) is called,
       but if you already have cache saved you can optionally pass it in
       to save time by not re-loading the datafiles.

 - swmf2vtk(filetag, use_ascii=False, cache=None)
    * writes a 3d_*.vtk for the inputed file tag 3d_* .
    * Note: If cache is None, then read_all(filetag) is called,
      but if you already have cache saved you can optionally pass it in
      to save time by not re-loading the datafiles.
