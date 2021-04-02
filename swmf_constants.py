#### Values from the SWMF soucecode, Batsrus.jl code, 
#### and also some extra named contants and indices


# Possible values for the status variable
Unset_     = -100 # index for unset values (that are otherwise larger)
Unused_      = -1 # unused block (not a leaf)
Refine_      = -2 # parent block to be refined
DontCoarsen  = -3 # block not to be coarsened
Coarsen_     = -4 # child block to be coarsened
Used_        =  1 # currently used block (leaf)
RefineNew_   =  2 # child block to be refined
Refined_     =  3 # refined child block
CoarsenNew_  =  4 # parent block to be coarsened
Coarsened_   =  5 # coarsened parent block

# Deepest AMR level relative to root nodes (limited by 32 bit integers)
MaxLevel = 30

'''from SWMF/GM/BATSRUS/srcBATL_/BATL_tree.f90
  ! The maximum integer coordinate for a given level below root nodes
  ! Implied do loop was not understooed by the pgf90 compiler, so list them
  integer, parameter, public :: MaxCoord_I(0:MaxLevel) = &
       [ 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, &
       16384, 32768, 65536, 131072, 262144, 524288, 1048576, 2097152, &
       4194304, 8388608, 16777216, 33554432, 67108864, 134217728, &
       268435456, 536870912, 1073741824 ]
'''#why are they indexing this one from 0?

MaxCoord_I =  [ 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192,
                16384, 32768, 65536, 131072, 262144, 524288, 1048576, 2097152, 
                4194304, 8388608, 16777216, 33554432, 67108864, 134217728,
                268435456, 536870912, 1073741824 ]

# check the copied hardcoded things are what I think they are
assert(len(MaxCoord_I) == MaxLevel+1)
assert(MaxCoord_I == [2**i for i in range(len(MaxCoord_I))])

# Named indexes of iTree_IA
Status_   =  1
Level_    =  2 # grid level
Proc_     =  3 # processor index
Block_    =  4 # block index
MinLevel_ =  5 # minimum level allowed
MaxLevel_ =  6 # maximum level allowed
Coord0_   =  6 # equal to Coord1_-1
Coord1_   =  7 # coordinate of node in 1st dimension
Coord2_   =  8 # coordinate of node in 2nd dimension
Coord3_   =  9 # coordinate of node in 3rd dimension
CoordLast_=  9 # Coord0_ + MaxDim (?)
Parent_   = 10 # Parent_ must be equal to Child0_
Child0_   = 10 #
Child1_   = Child0_ + 1
#ChildLast_= Child0_ + nChild

'''
the status of the node (used, unused, to be refined, to be coarsened, etc.);
the current, the maximum allowed and minimum allowed AMR levels for this node;
the three integer coordinates with respect to the whole grid;
the index of the parent node (if any);
the indexes of the children nodes (if any);
the processor index where the block is stored for active nodes;
the local block index for active nodes.
'''

# my added named constants
ROOTNODE_ = 1
