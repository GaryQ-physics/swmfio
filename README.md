# swmf_file_reader

reads the .out, .tree, and .info files for the data from the
(BATSRUS) magnetosphere component of an SWMF run.

note, currently uses spacepy to read the .out file,
and scipy.io.FortranFiles for the .tree file
but eventually plans to do so without a spacepy dependency.

Tested for python 2.7 and 3.7

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
