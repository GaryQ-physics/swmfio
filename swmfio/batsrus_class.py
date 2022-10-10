import numpy as np 
from numba import njit, types
from numba.typed import Dict
from numba.experimental import jitclass

import swmfio
from swmfio.constants import Used_,Unused_,Status_,Level_,Parent_,Child1_,Coord1_,CoordLast_
from swmfio.read_batsrus import read_tree, read_data
from swmfio.util import unravel_index

# npts -> nCell
# block_parent_id   -> node_parent_id   
# block_child_ids   -> node_child_ids   
# block_amr_levels  -> node_amr_levels  
# block_x_min       -> node_x_min       
# block_y_min       -> node_y_min       
# block_z_min       -> node_z_min       
# block_x_max       -> node_x_max       
# block_y_max       -> node_y_max       
# block_z_max       -> node_z_max       
# block_child_count -> node_child_count 

@njit
def F2P(fortran_index):
    return fortran_index - 1

@njit
def P2F(python_index):
    return python_index + 1

spec = [
            ('nDim'      , types.int32    ),
            ('nI'        , types.int32    ),
            ('nJ'        , types.int32    ),
            ('nK'        , types.int32    ),
            ('xGlobalMin', types.float32  ),
            ('yGlobalMin', types.float32  ),
            ('zGlobalMin', types.float32  ),
            ('xGlobalMax', types.float32  ),
            ('yGlobalMax', types.float32  ),
            ('zGlobalMax', types.float32  ),

            ('rootnode'          , types.int32      ),
            ('block_parent_id'   , types.int32[:]   ),
            ('block_child_ids'   , types.int32[:,:] ),
            ('block_amr_levels'  , types.int32[:]   ),
            ('block_x_min'       , types.float32[:] ),
            ('block_y_min'       , types.float32[:] ),
            ('block_z_min'       , types.float32[:] ),
            ('block_x_max'       , types.float32[:] ),
            ('block_y_max'       , types.float32[:] ),
            ('block_z_max'       , types.float32[:] ),
            ('block_child_count' , types.int8[:]    ),

            ('data_arr'  , types.float32[:,:]                              ),
            ('DataArray' , types.float32[:,:,:,:,:]                        ),
            ('varidx'    , types.DictType(types.unicode_type, types.int32) ),

            ('block2node', types.int32[:]     ),
            ('node2block', types.int32[:]     ),
            ('file'      , types.unicode_type )
        ]

#            ('iRatio_D'  , types.int32[:]                                 ),
#            ('nRoot_D'   , types.int32[:]                                 ),
#            ('nInfo'     , types.int32                                    ),

#@jitclass(spec)
class BatsrusClass:

    def __init__(self,
                    nDim      ,
                    nI        ,
                    nJ        ,
                    nK        ,
                    xGlobalMin,
                    yGlobalMin,
                    zGlobalMin,
                    xGlobalMax,
                    yGlobalMax,
                    zGlobalMax,

                    rootnode          ,
                    block_parent_id   ,
                    block_child_ids   ,
                    block_amr_levels  ,
                    block_x_min       ,
                    block_y_min       ,
                    block_z_min       ,
                    block_x_max       ,
                    block_y_max       ,
                    block_z_max       ,
                    block_child_count ,

                    data_arr     ,
                    DataArray    ,
                    varidx       ,

                    block2node   ,
                    node2block   ,
                    file):

        self.nDim              = nDim
        self.nI                = nI
        self.nJ                = nJ
        self.nK                = nK
        self.xGlobalMin        = xGlobalMin
        self.yGlobalMin        = yGlobalMin
        self.zGlobalMin        = zGlobalMin
        self.xGlobalMax        = xGlobalMax
        self.yGlobalMax        = yGlobalMax
        self.zGlobalMax        = zGlobalMax

        self.rootnode          = rootnode
        self.block_parent_id   = block_parent_id
        self.block_child_ids   = block_child_ids
        self.block_amr_levels  = block_amr_levels
        self.block_x_min       = block_x_min
        self.block_y_min       = block_y_min
        self.block_z_min       = block_z_min
        self.block_x_max       = block_x_max
        self.block_y_max       = block_y_max
        self.block_z_max       = block_z_max
        self.block_child_count = block_child_count 

        self.data_arr          = data_arr
        self.DataArray         = DataArray
        self.varidx            = varidx

        self.block2node        = block2node
        self.node2block        = node2block
        self.file              = file

        # map blocks <=> nodes for interpolation; add volume variable
        for iBlockP in range(block2node.size):
            #swmfio.logger.info(f"Working on block {iBlockP}/{block2node.size}")
            print(f"Working on block {iBlockP}/{block2node.size}")
            iNodeP = F2P( self.find_tree_node(data_arr[iBlockP*nI*nJ*nK, 0:3]) )
            print(f"Block {iBlockP} has iNodeP = {iNodeP}")
            self.block2node[iBlockP] = iNodeP
            self.node2block[iNodeP] = iBlockP

            epsilonX = DataArray[varidx['x'], 1,0,0,iBlockP] -  DataArray[varidx['x'], 0,0,0,iBlockP]
            epsilonY = DataArray[varidx['y'], 0,1,0,iBlockP] -  DataArray[varidx['y'], 0,0,0,iBlockP]
            epsilonZ = DataArray[varidx['z'], 0,0,1,iBlockP] -  DataArray[varidx['z'], 0,0,0,iBlockP]
            for k in range(nK):
                for j in range(nJ):
                    for i in range(nI):
                        self.DataArray[varidx['measure'], i, j, k, iBlockP] = epsilonX*epsilonY*epsilonZ


    def find_tree_node(self, point):

        xin = self.xGlobalMin <= point[0] <= self.xGlobalMax
        yin = self.yGlobalMin <= point[1] <= self.yGlobalMax
        zin = self.zGlobalMin <= point[2] <= self.zGlobalMax

        if not (xin and yin and zin): 
            raise RuntimeError('point out of simulation volume')

        iNode = self.rootnode

        n = 0;
        while True:
            if self.block_child_count[F2P(iNode)] == 0:
                print("Node does not have children. Done.")
                break

            print(f"n = {n}")
            n = n + 1
            for j in range(self.block_child_count[F2P(iNode)]):
                print(f"{j}/{self.block_child_count[F2P(iNode)]}")
                child = self.block_child_ids[j, F2P(iNode)]

                xin = self.block_x_min[F2P(child)] <= point[0] <= self.block_x_max[F2P(child)]
                yin = self.block_y_min[F2P(child)] <= point[1] <= self.block_y_max[F2P(child)]
                zin = self.block_z_min[F2P(child)] <= point[2] <= self.block_z_max[F2P(child)]

                if xin and yin and zin:
                    print(f"Found node for given point.")
                    print(f"x = {point[0]} block_x_min = {self.block_x_min[F2P(child)]} block_x_max = {self.block_x_max[F2P(child)]}")
                    print(f"y = {point[1]} block_y_min = {self.block_y_min[F2P(child)]} block_y_max = {self.block_y_max[F2P(child)]}")
                    print(f"z = {point[2]} block_z_min = {self.block_y_min[F2P(child)]} block_z_max = {self.block_y_max[F2P(child)]}")
                    iNode = child
                    break

        return iNode


    def interpolate(self, point, var):

        _x = self.varidx['x']
        _y = self.varidx['y']
        _z = self.varidx['z']
        iVar = self.varidx[var]

        DA = self.DataArray
        nVar, nI, nJ, nK, nBlock = DA.shape

        iNode = self.find_tree_node(point)
        iBlockP = self.node2block[F2P(iNode)]

        # get the gridspacing in x,y,z
        gridspacingX = DA[_x,1,0,0,iBlockP] - DA[_x,0,0,0,iBlockP]
        gridspacingY = DA[_y,0,1,0,iBlockP] - DA[_y,0,0,0,iBlockP]
        gridspacingZ = DA[_z,0,0,1,iBlockP] - DA[_z,0,0,0,iBlockP]

        # i0 is s.t. the highest index s.t. the x coordinate of the 
        #  corresponding cell block_data[iNode,i0,:,:]  is still less than point[0]
        i0 = (point[0] - DA[_x, 0, 0, 0, iBlockP])/gridspacingX
        j0 = (point[1] - DA[_y, 0, 0, 0, iBlockP])/gridspacingY
        k0 = (point[2] - DA[_z, 0, 0, 0, iBlockP])/gridspacingZ
        #if i0.is_integer() and j0.is_integer() and k0.is_integer(): # doesnt work in numba
        #    is_native = True
        #    print( (iBlockP,int(i0),int(j0),int(k0)) )# (iBlockP,i,j,k)
        #    i0 = int(i0)
        #    j0 = int(j0)
        #    k0 = int(k0)
        i0 = int(np.floor(i0))
        j0 = int(np.floor(j0))
        k0 = int(np.floor(k0))

        # i1 = i0+1 is the lowest index s.t. the x coordinate of the 
        #  corresponding cell block_data[iNode,i1,:,:]  is still greater than point[0]
        # together, i0 and i1 form the upper and lower bounds for a linear interpolation in x
        # likewise for j0,j1,y  and k0,k1,z

        # TODO: implement better interpolation at ends of block.
        # This method effectively makes it nearest neighbor at the ends
        if i0 == -1:
            i0 = 0
            i1 = 0
        elif i0 == nI-1:
            i1 = nI-1
        else:
            i1 = i0 + 1

        if j0 == -1:
            j0 = 0
            j1 = 0
        elif j0 == nJ-1:
            j1 = nJ-1
        else:
            j1 = j0 + 1

        if k0 == -1:
            k0 = 0
            k1 = 0
        elif k0 == nK-1:
            k1 = nK-1
        else:
            k1 = k0 + 1

        # all together i0,i1,j0, etc... form a cube of side length "gridpacing"
        # To do trilinear interpolation within, define xd as the distanc
        # along x of point within that cube, in units of "gridspacing"
        xd = (point[0] - DA[_x, i0, 0 , 0 , iBlockP])/gridspacingX
        yd = (point[1] - DA[_y, 0 , j0, 0 , iBlockP])/gridspacingY
        zd = (point[2] - DA[_z, 0 , 0 , k0, iBlockP])/gridspacingZ

        #https://en.wikipedia.org/wiki/Trilinear_interpolation
        c000 = DA[iVar,  i0, j0, k0,  iBlockP]
        c001 = DA[iVar,  i0, j0, k1,  iBlockP]
        c010 = DA[iVar,  i0, j1, k0,  iBlockP]
        c100 = DA[iVar,  i1, j0, k0,  iBlockP]
        c011 = DA[iVar,  i0, j1, k1,  iBlockP]
        c110 = DA[iVar,  i1, j1, k0,  iBlockP]
        c101 = DA[iVar,  i1, j0, k1,  iBlockP]
        c111 = DA[iVar,  i1, j1, k1,  iBlockP]

        c00 = c000*(1.-xd) + c100*xd
        c01 = c001*(1.-xd) + c101*xd
        c10 = c010*(1.-xd) + c110*xd
        c11 = c011*(1.-xd) + c111*xd

        c0 = c00*(1.-yd) + c10*yd
        c1 = c01*(1.-yd) + c11*yd

        c = c0*(1.-zd) + c1*zd
        return c


    def get_native_partial_derivatives(self, indx, var):
        _x = self.varidx['x']
        _y = self.varidx['y']
        _z = self.varidx['z']
        iVar = self.varidx[var]

        DA = self.DataArray
        nVar, nI, nJ, nK, nBlock = DA.shape

        i,j,k,iBlockP = unravel_index(indx, (nI,nJ,nK,nBlock), order='F')
        assert(indx == i + nI*j + nI*nJ*k + nI*nJ*nK*iBlockP) 

        partials = np.empty(3, dtype=np.float32)

        epsilonX = DA[_x,1,0,0,iBlockP] - DA[_x,0,0,0,iBlockP]
        epsilonY = DA[_y,0,1,0,iBlockP] - DA[_y,0,0,0,iBlockP]
        epsilonZ = DA[_z,0,0,1,iBlockP] - DA[_z,0,0,0,iBlockP]

        if i == 0:
            partials[0] = (DA[iVar, 1, j, k , iBlockP] - DA[iVar, 0, j, k, iBlockP])/(epsilonX)
        elif i == nI-1:
            partials[0] = (DA[iVar, nI-1, j, k, iBlockP] - DA[iVar, nI-2, j, k, iBlockP])/(epsilonX)
        else:
            partials[0] = (DA[iVar, i+1, j  , k, iBlockP] - DA[iVar, i-1, j, k, iBlockP])/(2*epsilonX)

        if j == 0:
            partials[1] = (DA[iVar, i, 1, k, iBlockP] - DA[iVar, i, 0, k, iBlockP])/(epsilonY)
        elif j == nJ-1:
            partials[1] = (DA[iVar, i, nJ-1, k, iBlockP] - DA[iVar, i, nJ-2, k, iBlockP])/(epsilonY)
        else:
            partials[1] = (DA[iVar, i, j+1, k, iBlockP] - DA[iVar, i, j-1, k, iBlockP])/(2*epsilonY)

        if k == 0:
            partials[2] = (DA[iVar, i, j, 1, iBlockP] - DA[iVar, i, j, 0, iBlockP])/(epsilonZ)
        elif k == nK-1:
            partials[2] = (DA[iVar, i, j, nK-1, iBlockP] - DA[iVar, i, j, nK-2, iBlockP])/(epsilonZ)
        else:
            partials[2] = (DA[iVar, i, j, k+1, iBlockP] - DA[iVar, i, j, k-1, iBlockP])/(2*epsilonZ)

        return partials


def get_class_from_native(file):

    iTree_IA, iRatio_D, nRoot_D, info = read_tree(file)
    data_arr, variables = read_data(file)

    nI = int(info['BlockSize1'])
    nJ = int(info['BlockSize2'])
    nK = int(info['BlockSize3'])
    xGlobalMin = float(info['Coord1Min'])
    yGlobalMin = float(info['Coord2Min'])
    zGlobalMin = float(info['Coord3Min'])
    xGlobalMax = float(info['Coord1Max'])
    yGlobalMax = float(info['Coord2Max'])
    zGlobalMax = float(info['Coord3Max'])

    nNode = iTree_IA.shape[1]
    npts = data_arr.shape[0]
    nBlock = npts//(nI*nJ*nK) if npts%(nI*nJ*nK)==0 else -1

    assert(not np.isfortran(data_arr))
    DataArray = data_arr.transpose()
    assert(np.isfortran(DataArray))
    DataArray = DataArray.reshape((len(variables)+1, nI, nJ, nK, nBlock), order='F')
    assert(np.isfortran(DataArray))

    varidx = Dict.empty(
        key_type=types.unicode_type,
        value_type=types.int32,
        )
    varidx['measure'] = np.int32(len(variables))
    for ivar,var in enumerate(variables):
        varidx[var] = np.int32(ivar)

    # In what follows, the P in iNodeP and iBlockP stands for Python-like
    # indexing (as opposed to Fortran)
    #
    # iNodeP indexes all nodes of the tree, from 0 to nNode-1,
    # and thus the "iNode" to be used in the other functions is simply iNodeP+1, or P2F(iNodeP)
    # 
    # iBlockP indexes all the blocks used, from 0 to nBlock-1. There is one for
    # each node with a status of used. 
    #
    # Note, nBlock*nI*nJ*nK = total number of batsrus cells (npts)

    # initialize arrays to -1 (invalid index), will be computed in __init__
    block2node = -np.ones((nBlock,), dtype=np.int32)
    node2block = -np.ones((nNode,), dtype=np.int32)

    # from SWMF/GM/BATSRUS/srcBATL/BATL_tree.f90 line 951 with substitutions
    def get_tree_position(iNode):
        '''
        Calculate normalized position of the edges of node iNode.
        Zero is at the minimum boundary of the grid, one is at the max boundary
        the tree is described by iTree_IA and pr
        '''
        iLevel = iTree_IA[F2P(Level_), F2P(iNode)]

        MaxIndex_D = ((2**(iLevel)-1)*(iRatio_D-1) + 1)*nRoot_D
        # Note: in the common case of iRatio_D=[2,2,2] and nRoot_D=[1,1,1]:
        # MaxIndex_D[all] = ((2**(iLevel)-1)*(2-1) + 1)*1
        #                 = 2**iLevel

        assert(MaxIndex_D.shape == (3,))
        assert(np.all( MaxIndex_D == 2**iLevel ))
        # note that if gridspacing = (256./8.)*0.5**iLevel, then gridspacing*MaxIndex_D[all] == 256./8. == 32. )

        # Convert to real by adding -1.0 or 0.0 for the two edges, respectively
        block_coords = iTree_IA[F2P(Coord1_):F2P(CoordLast_)+1,F2P(iNode)] # not seperatly defined in BATL_tree.f90
        PositionMin_D = (block_coords - 1.0)/MaxIndex_D
        PositionMax_D = (block_coords + 0.0)/MaxIndex_D

        return PositionMin_D, PositionMax_D


    block_parent_id = iTree_IA[Parent_, :].copy()
    block_child_ids = iTree_IA[F2P(Child1_):F2P(Child1_)+8, :].copy()
    block_amr_levels = iTree_IA[F2P(Level_), :].copy()

    block_x_min = np.empty(nNode, dtype=np.float32)
    block_y_min = np.empty(nNode, dtype=np.float32)
    block_z_min = np.empty(nNode, dtype=np.float32)
    block_x_max = np.empty(nNode, dtype=np.float32)
    block_y_max = np.empty(nNode, dtype=np.float32)
    block_z_max = np.empty(nNode, dtype=np.float32)
    block_child_count = np.empty(nNode, dtype=np.int8)
    for iNodeP in range(nNode):
        if   iTree_IA[F2P(Status_), iNodeP] == Used_:
            block_child_count[iNodeP] = 0
        elif iTree_IA[F2P(Status_), iNodeP] == Unused_:
            block_child_count[iNodeP] = 8
        else:
            assert(False)

        x_start = xGlobalMin
        y_start = yGlobalMin
        z_start = zGlobalMin
        x_range = xGlobalMax - xGlobalMin
        y_range = yGlobalMax - yGlobalMin
        z_range = zGlobalMax - zGlobalMin

        iLevel = iTree_IA[F2P(Level_), iNodeP]
        assert(nI == nJ == nK) #!!!!
        assert(x_range == y_range == z_range)
        gridspacing = (x_range/nI)*0.5**iLevel

        PositionMin_D, PositionMax_D = get_tree_position(P2F(iNodeP))
        block_x_min[iNodeP] = x_range*(PositionMin_D[0]) + x_start
        block_y_min[iNodeP] = y_range*(PositionMin_D[1]) + y_start
        block_z_min[iNodeP] = z_range*(PositionMin_D[2]) + z_start
        block_x_max[iNodeP] = x_range*(PositionMax_D[0]) + x_start
        block_y_max[iNodeP] = y_range*(PositionMax_D[1]) + y_start
        block_z_max[iNodeP] = z_range*(PositionMax_D[2]) + z_start


    batsclass = BatsrusClass(
                      nDim              = 3         ,
                      nI                = nI        ,
                      nJ                = nJ        ,
                      nK                = nK        ,
                      xGlobalMin        = xGlobalMin,
                      yGlobalMin        = yGlobalMin,
                      zGlobalMin        = zGlobalMin,
                      xGlobalMax        = xGlobalMax,
                      yGlobalMax        = yGlobalMax,
                      zGlobalMax        = zGlobalMax,

                      rootnode          = 1,
                      block_parent_id   = block_parent_id   ,
                      block_child_ids   = block_child_ids   ,
                      block_amr_levels  = block_amr_levels  ,
                      block_x_min       = block_x_min       ,
                      block_y_min       = block_y_min       ,
                      block_z_min       = block_z_min       ,
                      block_x_max       = block_x_max       ,
                      block_y_max       = block_y_max       ,
                      block_z_max       = block_z_max       ,
                      block_child_count = block_child_count ,

                      data_arr          = data_arr     ,
                      DataArray         = DataArray    ,
                      varidx            = varidx       ,

                      block2node        = block2node   ,
                      node2block        = node2block   ,
                      file              = file
                )

    return batsclass


def get_class_from_cdf(file):

    import cdflib.cdfread as cdfread

    cdf = cdfread.CDF(file)
    globatts = cdf.globalattsget()

    npts = int(globatts['number_of_cells'])
    nBlock = int(globatts['number_of_blocks'])
    nI = int(globatts['special_parameter_NX'])
    nJ = int(globatts['special_parameter_NY'])
    nK = int(globatts['special_parameter_NZ'])
    assert( nBlock*nI*nJ*nK == npts )

    swmfio.logger.info(f"npts = {npts}")
    swmfio.logger.info(f"nBlock = {nBlock}")
    swmfio.logger.info(f"nI/nJ/nK = {nI}/{nJ}/{nK}")

    nNode = cdf.varget('block_amr_levels').size
    block2node = -np.ones((nBlock,), dtype=np.int32)
    node2block = -np.ones((nNode,), dtype=np.int32)

    swmfio.logger.info(f"nNode = {nNode}")

    varidx = Dict.empty(key_type=types.unicode_type, value_type=types.int64,)
    units = {}

    nVar = 0
    for cdfvar in cdf.cdf_info()['zVariables']:
        if cdf.varget(cdfvar).shape == (1, npts):
            nVar += 1

    nVar += 1 # for added measure (volume) variable

    data_arr = np.empty((npts,nVar), dtype=np.float32);
    data_arr[:,:] = np.nan

    iVar = 0
    for cdfvar in cdf.cdf_info()['zVariables']:
        try:
            var = cdf.varattsget(cdfvar)['Original Name']
        except:
            var = cdfvar

        if cdf.varget(cdfvar).shape == (1, npts):
            #swmfio.logger.info(f"Reading = {cdfvar}")
            data_arr[:, iVar] = cdf.varget(cdfvar)[0,:]
            units[var] = cdf.varattsget(cdfvar)['units']
            varidx[var] = iVar
            iVar += 1

    swmfio.logger.info("varidx = {}".format(varidx))

    varidx['measure'] = iVar

    assert(not np.isfortran(data_arr))

    DataArray = data_arr.transpose()
    assert(np.isfortran(DataArray))

    DataArray = DataArray.reshape((nVar, nI, nJ, nK, nBlock), order='F')
    assert(np.isfortran(DataArray))

    block_child_ids = np.array([
                                    cdf.varget('block_child_id_1')[0,:],
                                    cdf.varget('block_child_id_2')[0,:],
                                    cdf.varget('block_child_id_3')[0,:],
                                    cdf.varget('block_child_id_4')[0,:],
                                    cdf.varget('block_child_id_5')[0,:],
                                    cdf.varget('block_child_id_6')[0,:],
                                    cdf.varget('block_child_id_7')[0,:],
                                    cdf.varget('block_child_id_8')[0,:],
                                ])
    block_child_ids = np.array(P2F(block_child_ids), dtype=np.int32)

    batsclass = BatsrusClass(
                      nDim              = globatts['grid_system_1_number_of_dimensions'],
                      nI                = nI,
                      nJ                = nJ,
                      nK                = nK,
                      xGlobalMin        = globatts['global_x_min'],
                      yGlobalMin        = globatts['global_y_min'],
                      zGlobalMin        = globatts['global_z_min'],
                      xGlobalMax        = globatts['global_x_max'],
                      yGlobalMax        = globatts['global_y_max'],
                      zGlobalMax        = globatts['global_z_max'],

                      rootnode          = P2F( cdf.varget('block_at_amr_level')[0,0] ),
                      block_parent_id   = cdf.varget('block_parent_id')[0,:],
                      block_child_ids   = block_child_ids,
                      block_amr_levels  = np.array(cdf.varget('block_amr_levels')[0,:], dtype=np.int32),
                      block_x_min       = cdf.varget('block_x_min')[0,:],
                      block_y_min       = cdf.varget('block_y_min')[0,:],
                      block_z_min       = cdf.varget('block_z_min')[0,:],
                      block_x_max       = cdf.varget('block_x_max')[0,:],
                      block_y_max       = cdf.varget('block_y_max')[0,:],
                      block_z_max       = cdf.varget('block_z_max')[0,:],
                      block_child_count = np.array(cdf.varget('block_child_count')[0,:], dtype=np.int8),

                      data_arr          = data_arr     ,
                      DataArray         = DataArray    ,
                      varidx            = varidx       ,

                      block2node        = block2node   ,
                      node2block        = node2block   ,
                      file              = file
                )

    return batsclass
