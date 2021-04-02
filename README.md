# swmf_file_reader

reads the .out, .tree, and .info files for the data from the
(BATSRUS) magnetosphere component of an SWMF run.

note, currently uses spacepy to read the .out file,
and scipy.io.FortranFiles for the .tree file
but eventually plans to do so without a spacepy dependency.

note, currently writing vtk only works with python 2
