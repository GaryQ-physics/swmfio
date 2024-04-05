
The batsrus_interoplator uses the CCMC OCTREE library. From Kamodo, download 
the OCTREE library:

https://github.com/nasa/Kamodo/tree/master/kamodo_ccmc/readers/OCTREE_BLOCK_GRID

The OCTREE files must be located with the other swmfio Python files.  So create
{Python Code}/swmfio/swmfio/OCTREE_BLOCK_GRID, where {Python Code}/swmfio/swmfio 
contains the swmfio Python files.  The directory {Python Code}/swmfio/swmfio 
should contain:

BATSRUS_interpolator.py	
batsrus_imperative.py	
util.py
OCTREE_BLOCK_GRID	
constants.py		
vtk_export.py
__init__.py		
read_batsrus.py		
write_vtk.py
__pycache__		
read_rim.py
batsrus_class.py	test

To build the CCMC OCTREE library, from the terminal, change into the 
OCTREE_BLOCK_GRID directory, and execute:

python interpolate_amrdata_extension_build.py

To demo the interpolator, execute demo_interpolator.py




