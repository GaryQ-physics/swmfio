import numpy as np
from swmf_file_reader.swmf_constants import Used_,Status_,Level_,Parent_,Child0_,Child1_,Coord1_,CoordLast_,ROOTNODE_
from swmf_file_reader.read_batsrus import read_tree_file, read_out_file
from swmf_file_reader import util

from numba import njit, types
from numba.typed import Dict
from numba.experimental import jitclass

@njit
def F2P(fortran_index):
    return fortran_index - 1

@njit
def P2F(python_index):
    return python_index + 1

spec = [
            ('nDim'      , types.int32                                    ),
            ('nI'        , types.int32                                    ),
            ('nJ'        , types.int32                                    ),
            ('nK'        , types.int32                                    ),
            ('xGlobalMin', types.float32                                  ),
            ('yGlobalMin', types.float32                                  ),
            ('zGlobalMin', types.float32                                  ),
            ('xGlobalMax', types.float32                                  ),
            ('yGlobalMax', types.float32                                  ),
            ('zGlobalMax', types.float32                                  ),
            ('nInfo'     , types.int32                                    ),
            ('nNode'     , types.int32                                    ),
            ('iRatio_D'  , types.int32[:]                                 ),
            ('nRoot_D'   , types.int32[:]                                 ),

            ('data_arr'  , types.float32[:,:]                             ),
            ('DataArray' , types.float32[:,:,:,:,:]                       ),
            ('iTree_IA'  , types.int32[:,:]                               ),
            ('block2node', types.int32[:]                               ),
            ('node2block', types.int32[:]                               ),
            ('varidx'    , types.DictType(types.unicode_type, types.int32)),
            ]
@jitclass(spec)
class BatsrusClass:
    def __init__(self, 
                    nDim       ,
                    nI         ,
                    nJ         ,
                    nK         ,
                    xGlobalMin ,
                    yGlobalMin ,
                    zGlobalMin ,
                    xGlobalMax ,
                    yGlobalMax ,
                    zGlobalMax ,
                    nInfo      ,
                    nNode      ,
                    iRatio_D   ,
                    nRoot_D    ,

                    data_arr   ,
                    DataArray  ,
                    iTree_IA   ,
                    block2node ,
                    node2block ,
                    varidx         ):

        self.nDim       = nDim      
        self.nI         = nI        
        self.nJ         = nJ        
        self.nK         = nK        
        self.xGlobalMin = xGlobalMin
        self.yGlobalMin = yGlobalMin
        self.zGlobalMin = zGlobalMin
        self.xGlobalMax = xGlobalMax
        self.yGlobalMax = yGlobalMax
        self.zGlobalMax = zGlobalMax
        self.nInfo      = nInfo     
        self.nNode      = nNode     
        self.iRatio_D   = iRatio_D  
        self.nRoot_D    = nRoot_D   

        self.data_arr   = data_arr     
        self.DataArray  = DataArray    
        self.iTree_IA   = iTree_IA     
        self.block2node = block2node   
        self.node2block = node2block   
        self.varidx     = varidx       

        npts = data_arr.shape[0]
        nBlock = npts//(nI*nJ*nK) if npts%(nI*nJ*nK)==0 else -1
        # map blocks and nodes
        for iBlockP in range(nBlock):
            iNodeP = F2P( self.find_tree_node(data_arr[iBlockP*nI*nJ*nK, 0:3]) )
            self.block2node[iBlockP] = iNodeP
            self.node2block[iNodeP] = iBlockP

        for iBlockP in range(nBlock):
            epsilonX = DataArray[varidx['x'], 1,0,0,iBlockP] -  DataArray[varidx['x'], 0,0,0,iBlockP]
            epsilonY = DataArray[varidx['y'], 0,1,0,iBlockP] -  DataArray[varidx['y'], 0,0,0,iBlockP]
            epsilonZ = DataArray[varidx['z'], 0,0,1,iBlockP] -  DataArray[varidx['z'], 0,0,0,iBlockP]
            for k in range(nK):
                for j in range(nJ):
                    for i in range(nI):
                        self.DataArray[varidx['measure'], i,j,k,iBlockP] = epsilonX*epsilonY*epsilonZ

    # from SWMF/GM/BATSRUS/srcBATL/BATL_tree.f90 line 951 with substitutions
    def get_tree_position(self, iNode):
        '''
        Calculate normalized position of the edges of node iNode.
        Zero is at the minimum boundary of the grid, one is at the max boundary
        the tree is described by iTree_IA and pr
        '''
        iLevel = self.iTree_IA[F2P(Level_), F2P(iNode)]

        MaxIndex_D = ((2**(iLevel)-1)*(self.iRatio_D-1) + 1)*self.nRoot_D
        # Note: in the common case of iRatio_D=[2,2,2] and nRoot_D=[1,1,1]:
        # MaxIndex_D[all] = ((2**(iLevel)-1)*(2-1) + 1)*1
        #                 = 2**iLevel

        assert(MaxIndex_D.shape == (3,))
        assert(np.all( MaxIndex_D == 2**iLevel ))
        # note that if gridspacing = (256./8.)*0.5**iLevel, then gridspacing*MaxIndex_D[all] == 256./8. == 32. )

        # Convert to real by adding -1.0 or 0.0 for the two edges, respectively
        block_coords = self.iTree_IA[F2P(Coord1_):F2P(CoordLast_)+1,F2P(iNode)] # not seperatly defined in BATL_tree.f90
        PositionMin_D = (block_coords - 1.0)/MaxIndex_D
        PositionMax_D = (block_coords + 0.0)/MaxIndex_D

        return PositionMin_D, PositionMax_D

    def get_physical_dimensions(self, iNode, returnCenters=False):
        x_start = self.xGlobalMin
        y_start = self.yGlobalMin
        z_start = self.zGlobalMin
        x_range = self.xGlobalMax - self.xGlobalMin
        y_range = self.yGlobalMax - self.yGlobalMin
        z_range = self.zGlobalMax - self.zGlobalMin

        iLevel = self.iTree_IA[F2P(Level_), F2P(iNode)]
        assert(self.nI == self.nJ == self.nK)#!!!!
        assert(x_range == y_range == z_range)
        gridspacing = (x_range/self.nI)*0.5**iLevel

        PositionMin_D, PositionMax_D = self.get_tree_position(iNode)
        xmin = x_range*(PositionMin_D[0]) + x_start
        ymin = y_range*(PositionMin_D[1]) + y_start
        zmin = z_range*(PositionMin_D[2]) + z_start
        xmax = x_range*(PositionMax_D[0]) + x_start
        ymax = y_range*(PositionMax_D[1]) + y_start
        zmax = z_range*(PositionMax_D[2]) + z_start

        if returnCenters:
            xlims = (xmin+gridspacing/2., xmax-gridspacing/2.) #!! faster to use arrays??
            ylims = (ymin+gridspacing/2., ymax-gridspacing/2.)
            zlims = (zmin+gridspacing/2., zmax-gridspacing/2.)
            return xlims, ylims, zlims, gridspacing
        else:
            xminmax = (xmin, xmax)
            yminmax = (ymin, ymax)
            zminmax = (zmin, zmax)
            return xminmax, yminmax, zminmax, gridspacing

    # supposed to reproduce SWMF/GM/BATSRUS/srcBATL/BATL_tree.f90 line 975, but differently
    def find_tree_node(self, point):
        xin = self.xGlobalMin <= point[0] <= self.xGlobalMax
        yin = self.yGlobalMin <= point[1] <= self.yGlobalMax
        zin = self.zGlobalMin <= point[2] <= self.zGlobalMax

        if not (xin and yin and zin): 
            raise RuntimeError ('point out of simulation volume')

        iNode = ROOTNODE_
        while(True):
            if Used_ == self.iTree_IA[F2P(Status_), F2P(iNode)]:
                break

            for j in range(8):
                child = self.iTree_IA[F2P(Child1_+j), F2P(iNode)]
                xminmax, yminmax, zminmax, gridspacing = self.get_physical_dimensions(child,
                                                                                returnCenters=False)

                xin = xminmax[0] <= point[0] <= xminmax[1]
                yin = yminmax[0] <= point[1] <= yminmax[1]
                zin = zminmax[0] <= point[2] <= zminmax[1]

                if xin and yin and zin: 
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

        #TODO implement better interpolation at ends of block.
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

        # all together i0,i1,j0,ect... form a cube of side length "gridpacing" to do trililear interpolation within
        # define xd as the distance along x of point within that cube, in units of "gridspacing"
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

        i,j,k,iBlockP = util.unravel_index(indx, (nI,nJ,nK,nBlock), order='F')
        print(indx == i + nI*j + nI*nJ*k + nI*nJ*nK*iBlockP)

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


def return_class(filetag):
    iTree_IA, pr = read_tree_file(filetag)
    data_arr, variables = read_out_file(filetag)

    ## in what follows:
    #  the P in iNodeP and iBlockP stands for python like indexing (as oposed to fortran)
    #  
    #  iNodeP indexes all nodes of the tree, from 0 to nNode-1,
    #  and thus the "iNode" to be used in the other functions is simply iNodeP+1, or P2F(iNodeP)
    # 
    #  iBlockP indexes all the blocks used, from 0 to nBlock-1.
    #  There is one for each node with a status of used. 
    #  Note, nBlock*nI*nJ*nK = total number of batsrus cells (npts)
    npts = data_arr.shape[0]
    nI, nJ, nK = pr.nI, pr.nJ, pr.nK
    nBlock = npts//(nI*nJ*nK) if npts%(nI*nJ*nK)==0 else -1
    # initialize arrays to -1 (invalid index), will be computed in __init__
    block2node = -np.ones((nBlock,), dtype=np.int32)
    node2block = -np.ones((pr.nNode,), dtype=np.int32)

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

    batsclass = BatsrusClass(
                      nDim       = pr.nDim      ,
                      nI         = pr.nI        ,
                      nJ         = pr.nJ        ,
                      nK         = pr.nK        ,
                      xGlobalMin = pr.xGlobalMin,
                      yGlobalMin = pr.yGlobalMin,
                      zGlobalMin = pr.zGlobalMin,
                      xGlobalMax = pr.xGlobalMax,
                      yGlobalMax = pr.yGlobalMax,
                      zGlobalMax = pr.zGlobalMax,
                      nInfo      = pr.nInfo     ,
                      nNode      = pr.nNode     ,
                      iRatio_D   = pr.iRatio_D  ,
                      nRoot_D    = pr.nRoot_D   ,

                      data_arr   = data_arr     ,
                      DataArray  = DataArray    ,
                      iTree_IA   = iTree_IA     ,
                      block2node = block2node   ,
                      node2block = node2block   ,
                      varidx     = varidx       ,
                      )

    return batsclass


def main():
    cls = return_class('/tmp/3d__var_2_e20190902-041000-000')
    print(cls.varidx)
    print(cls.data_arr)

    onpoint = np.array([-146.,  -14.,  -14.])
    offpoint = np.array([-143.,  -15.,  -15.])
    print(cls.interpolate(onpoint, 'rho'))
    print(cls.interpolate(onpoint, 'x'))
    print(cls.interpolate(onpoint, 'y'))
    print(cls.interpolate(onpoint, 'z'))
    print(cls.interpolate(offpoint, 'rho'))
    print(cls.interpolate(offpoint, 'x'))
    print(cls.interpolate(offpoint, 'y'))
    print(cls.interpolate(offpoint, 'z'))
    print(cls.get_native_partial_derivatives(123456, 'rho'))
    print(cls.get_native_partial_derivatives(125356, 'rho'))
    print(cls.get_native_partial_derivatives(143456, 'rho'))
    print(cls.get_native_partial_derivatives(143456, 'x'))
    print(cls.get_native_partial_derivatives(143456, 'y'))
    print(cls.get_native_partial_derivatives(143456, 'z'))

if __name__ == '__main__':
    main()
