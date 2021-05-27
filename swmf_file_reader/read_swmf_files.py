from collections import namedtuple
import numpy as np
import scipy.io as sio
from swmf_file_reader.swmf_constants import Used_,Status_,Level_,Parent_,Child0_,Child1_,Coord1_,CoordLast_,ROOTNODE_
from swmf_file_reader.vtk_export_copy import vtk_export
from swmf_file_reader.named_var_indexes import nVarNeeded, index2str, str2index,str2index_typed, _x,_y,_z,_bx,_by,_bz

try:
    from numba import njit
except:
    print("\n\nWARNING: numba couldn't be imported so this will run slower\n\n")
    # construct a trivial decorator for njit 
    def njit(foo):
        return foo

BatsProps = namedtuple('BATSRUSProperties',['nDim','nI','nJ','nK', 
    'xGlobalMin','yGlobalMin','zGlobalMin','xGlobalMax','yGlobalMax','zGlobalMax',
    'nInfo','nNode','iRatio_D','nRoot_D'])

@njit
def F2P(fortran_index):
    return fortran_index - 1

@njit
def P2F(python_index):
    return python_index + 1


def read_info_file(filetag):
    info = {'filetag' : filetag}
    with open(filetag+'.info','r') as f:
        for line in f.readlines():
            if line == '\n' : continue
            if line[0] == '#': continue
            splt = line.split()
            if len(splt) == 2:
                info[splt[1]] = splt[0]
    return info


def read_tree_file(filetag):
    # first read info file
    info = read_info_file(filetag)

    ## load tree file
    ff = sio.FortranFile(filetag+".tree", 'r')
    nDim, nInfo, nNode = ff.read_ints(dtype=np.int32)
    iRatio_D = ff.read_ints(dtype=np.int32) # Array of refinement ratios
    nRoot_D = ff.read_ints(dtype=np.int32) # The number of root nodes in all dimension
    iTree_IA = ff.read_ints(dtype=np.int32).reshape((nInfo,nNode), order='F')

    pr = BatsProps( nInfo    = nInfo                      ,
                    nNode    = nNode                      ,
                    iRatio_D = iRatio_D                   ,
                    nRoot_D  = nRoot_D                    ,
                    nDim = int(info['nDim'])              ,
                    nI = int(info['BlockSize1'])          ,
                    nJ = int(info['BlockSize2'])          ,
                    nK = int(info['BlockSize3'])          ,
                    xGlobalMin = float(info['Coord1Min']) ,
                    yGlobalMin = float(info['Coord2Min']) ,
                    zGlobalMin = float(info['Coord3Min']) ,
                    xGlobalMax = float(info['Coord1Max']) ,
                    yGlobalMax = float(info['Coord2Max']) ,
                    zGlobalMax = float(info['Coord3Max'])       )

    ########################### check_thing_work #######################
    assert(pr.nDim == nDim)
    assert(iTree_IA.shape[1] == nNode)
    # Maximum number of ghost cells set by Config.pl script.
    # Valid values are 0,1,2,3,4,5
    nG = 2
    # Refinement ratios in the 3 dimensions. Either 1 or 2.
    # The values are set by the Config.pl script.
    iRatio, jRatio, kRatio = min(2, pr.nI), min(2, pr.nJ), min(2, pr.nK)
    # Number of dimensions in which grid adaptation is done
    nDimAmr = iRatio + jRatio + kRatio - 3
    assert(nDimAmr == nDim)
    assert(nDim == 3)
    # Number of children per node
    nChild = 2**nDimAmr
    assert(np.isfortran(iTree_IA))
    assert(np.all(iRatio_D == np.array([iRatio, jRatio, kRatio])))
    ####################################################################

    return iTree_IA, pr


def read_out_file(filetag):
    nDim = 3
    ff = sio.FortranFile(filetag+".out", 'r')

    header = ff.read_ints(dtype=np.uint8).tobytes().decode('UTF-8')
    nStep, Time, nDimOut, nParam, nVar = ff.read_ints(dtype=np.int32)
    n_D = ff.read_ints(dtype=np.int32)
    whatISthis=ff.read_ints(dtype=np.uint8).tobytes()
    try:
        print(whatISthis.decode('UTF-8'))
    except:
        pass
    allvariables = ff.read_ints(dtype=np.uint8).tobytes().decode('UTF-8')
    allvariables = allvariables.strip().split(' ')
    npts = n_D[0]; assert(n_D[1]==1 and n_D[2]==1 and n_D.size==3)

    #x,y,and z  are stored as one larger array in file. All others are stored seperately
    #as a results xyz count only once towards nVar.
    data_arr = np.empty((nVar+2, npts), dtype=np.float32); data_arr[:,:]=np.nan; assert(nVar+2==nVarNeeded)
    data_arr[0:3, :] = ff.read_reals(dtype=np.float32).reshape(3, npts) 
    for i in range(3, nVar+2):
        A = ff.read_reals(dtype=np.float32)
        data_arr[i, :] = A

    status = ff.read_reals(dtype=np.float32) #? whatever this means

    # check to make sure ff is at the end of the file
    #end = ff.read_reals(dtype=np.float32)
    #print(end)
    #exit()
    try:
        ff.read_reals(dtype=np.float32)
        assert(False)
    except TypeError: # nessesary for python 2 compatibility, since it seems that there is no sio._fortran.FortranEOFError without it
        pass
    except sio._fortran.FortranEOFError:
        pass

    expectedheader = "R R R Mp/cc km/s km/s km/s J/m3 nT nT nT nT nT nT nPa uA/m2 uA/m2 uA/m2 --"
    assert(expectedheader == header.strip())

    return data_arr, tuple(allvariables[:nVar+2])

# from SWMF/GM/BATSRUS/srcBATL/BATL_tree.f90 line 951 with substitutions
@njit
def get_tree_position(iNode, iTree_IA, pr):
    '''
    Calculate normalized position of the edges of node iNode.
    Zero is at the minimum boundary of the grid, one is at the max boundary
    the tree is described by iTree_IA and pr
    '''
    iLevel = iTree_IA[F2P(Level_), F2P(iNode)]

    MaxIndex_D = ((2**(iLevel)-1)*(pr.iRatio_D-1) + 1)*pr.nRoot_D
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

@njit
def get_physical_dimensions(iNode, iTree_IA, pr, returnCenters=False):
    x_start = pr.xGlobalMin
    y_start = pr.yGlobalMin
    z_start = pr.zGlobalMin
    x_range = pr.xGlobalMax - pr.xGlobalMin
    y_range = pr.yGlobalMax - pr.yGlobalMin
    z_range = pr.zGlobalMax - pr.zGlobalMin

    iLevel = iTree_IA[F2P(Level_), F2P(iNode)]
    assert(pr.nI == pr.nJ == pr.nK)
    assert(x_range == y_range == z_range)
    gridspacing = (x_range/pr.nI)*0.5**iLevel

    PositionMin_D, PositionMax_D = get_tree_position(iNode, iTree_IA, pr)
    xmin = x_range*(PositionMin_D[0]) + x_start
    ymin = y_range*(PositionMin_D[1]) + y_start
    zmin = z_range*(PositionMin_D[2]) + z_start
    xmax = x_range*(PositionMax_D[0]) + x_start
    ymax = y_range*(PositionMax_D[1]) + y_start
    zmax = z_range*(PositionMax_D[2]) + z_start

    if returnCenters:
        xlims = (xmin+gridspacing/2., xmax-gridspacing/2.)
        ylims = (ymin+gridspacing/2., ymax-gridspacing/2.)
        zlims = (zmin+gridspacing/2., zmax-gridspacing/2.)
        return xlims, ylims, zlims, gridspacing
    else:
        xminmax = (xmin, xmax)
        yminmax = (ymin, ymax)
        zminmax = (zmin, zmax)
        return xminmax, yminmax, zminmax, gridspacing

# supposed to reproduce SWMF/GM/BATSRUS/srcBATL/BATL_tree.f90 line 975, but differently
@njit
def find_tree_node(point, iTree_IA, pr):
    xin = pr.xGlobalMin <= point[0] <= pr.xGlobalMax
    yin = pr.yGlobalMin <= point[1] <= pr.yGlobalMax
    zin = pr.zGlobalMin <= point[2] <= pr.zGlobalMax

    if not (xin and yin and zin): 
        raise RuntimeError ('point out of simulation volume')

    iNode = ROOTNODE_
    while(True):
        if Used_ == iTree_IA[F2P(Status_), F2P(iNode)]:
            break

        for j in range(8):
            child = iTree_IA[F2P(Child1_+j), F2P(iNode)]
            xminmax, yminmax, zminmax, gridspacing = get_physical_dimensions(child, iTree_IA, pr,
                                                                            returnCenters=False)

            xin = xminmax[0] <= point[0] <= xminmax[1]
            yin = yminmax[0] <= point[1] <= yminmax[1]
            zin = zminmax[0] <= point[2] <= zminmax[1]

            if xin and yin and zin: 
                iNode = child
                break

    return iNode

@njit
def map_nodes_to_data(data_arr, iTree_IA, pr):
    # in what follows:
    #  the P in iNodeP and iBlockP stands for python like indexing (as oposed to fortran)
    #  
    #  iNodeP indexes all nodes of the tree, from 0 to nNode-1,
    #  and thus the "iNode" to be used in the other functions is simply iNodeP+1, or P2F(iNodeP)
    # 
    #  iBlockP indexes all the blocks used, from 0 to nBlock-1.
    #  There is one for each node with a status of used. 
    #  Note, nBlock*nI*nJ*nK = total number of batsrus cells (npts)
    npts = data_arr.shape[1]
    nI, nJ, nK = pr.nI, pr.nJ, pr.nK
    nBlock = npts//(nI*nJ*nK) if npts%(nI*nJ*nK)==0 else -1

    block2node = -np.ones((nBlock,), dtype=np.int64)
    node2block = -np.ones((pr.nNode,), dtype=np.int64)
    for iBlockP in range(nBlock):
        iNodeP = F2P( find_tree_node(data_arr[0:3, iBlockP*8**3], iTree_IA, pr) )
        block2node[iBlockP] = iNodeP
        node2block[iNodeP] = iBlockP
    return block2node, node2block


def read_all(filetag):
    iTree_IA, pr = read_tree_file(filetag)
    data_arr, variables = read_out_file(filetag)
    block2node, node2block = map_nodes_to_data(data_arr, iTree_IA, pr)

    assert(not np.isfortran(data_arr))
    DataArray = data_arr.reshape((data_arr.shape[0], block2node.size, pr.nK, pr.nJ, pr.nI))
    assert(not np.isfortran(DataArray))
    DataArray = DataArray.transpose((0,1,4,3,2)) # probably duplicates in memory, but without we'd need by (nVar,nBlock, nK,nJ,nI) 
    assert(not np.isfortran(DataArray))
    # would be better to have DataArray Fortran ordered with shape (nI,nJ,nK,nBlock,nVar), cause then we wouldnt need to rewrite
    # this isnt compatible at the moment though

    cache = {}
    cache['DataArray'] = DataArray
    cache['iTree_IA'] = iTree_IA
    cache['block2node'] = block2node
    cache['node2block'] = node2block
    cache['nBlock'] = DataArray.shape[1]
    cache['nI'] = DataArray.shape[2]
    cache['nJ'] = DataArray.shape[3]
    cache['nK'] = DataArray.shape[4]
    cache['filetag'] = filetag
    cache['pr'] = pr
    return cache

def find_index(filetag, point, cache=None, debug=False): # depreciated
    print('find_index DEPRECIATED')
    if cache is None:
        cache = read_all(filetag)
    else:
        assert(cache['filetag'] == filetag)

    def getvar(_var, iNode, i,j,k):
        return cache['DataArray'][str2index[_var],:,:,:,:][cache['node2block'][F2P(iNode)],i,j,k]

    iNode = find_tree_node(point, cache['iTree_IA'], cache['pr'])

    # get the gridspacing in x,y,z
    gridspacingX = getvar('x', iNode, 1,0,0) - getvar('x', iNode, 0,0,0)
    gridspacingY = getvar('y', iNode, 0,1,0) - getvar('y', iNode, 0,0,0)
    gridspacingZ = getvar('z', iNode, 0,0,1) - getvar('z', iNode, 0,0,0)

    # i0 is s.t. the highest index s.t. the x coordinate of the 
    #  corresponding cell block_data[iNode,i0,:,:]  is still less than point[0]
    i0 = (point[0] - getvar('x', iNode, 0, 0, 0))/gridspacingX
    j0 = (point[1] - getvar('y', iNode, 0, 0, 0))/gridspacingY
    k0 = (point[2] - getvar('z', iNode, 0, 0, 0))/gridspacingZ
    if i0.is_integer() and j0.is_integer() and k0.is_integer():
        return (cache['node2block'][F2P(iNode)],int(i0),int(j0),int(k0)) # (iBlockP,i,j,k)
    else:
        return None


@njit
def interpolate_jit(point, iVar, DA, iTree_IA, node2block, pr):
    _x,_y,_z = 0,1,2
    nI, nJ, nK = pr.nI, pr.nJ, pr.nK

    iNode = find_tree_node(point, iTree_IA, pr)
    iBlockP = node2block[F2P(iNode)]

    gridspacingX = DA[_x,iBlockP, 1,0,0] - DA[_x,iBlockP, 0,0,0]
    gridspacingY = DA[_y,iBlockP, 0,1,0] - DA[_y,iBlockP, 0,0,0]
    gridspacingZ = DA[_z,iBlockP, 0,0,1] - DA[_z,iBlockP, 0,0,0]

    # i0 is s.t. the highest index s.t. the x coordinate of the 
    #  corresponding cell block_data[iNode,i0,:,:]  is still less than point[0]
    i0 = (point[0] - DA[_x,iBlockP, 0, 0, 0])/gridspacingX
    j0 = (point[1] - DA[_y,iBlockP, 0, 0, 0])/gridspacingY
    k0 = (point[2] - DA[_z,iBlockP, 0, 0, 0])/gridspacingZ
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
    xd = (point[0] - DA[_x,iBlockP, i0, 0 , 0 ])/gridspacingX
    yd = (point[1] - DA[_y,iBlockP, 0 , j0, 0 ])/gridspacingY
    zd = (point[2] - DA[_z,iBlockP, 0 , 0 , k0])/gridspacingZ

    #https://en.wikipedia.org/wiki/Trilinear_interpolation
    c000 = DA[iVar,iBlockP,  i0, j0, k0]
    c001 = DA[iVar,iBlockP,  i0, j0, k1]
    c010 = DA[iVar,iBlockP,  i0, j1, k0]
    c100 = DA[iVar,iBlockP,  i1, j0, k0]
    c011 = DA[iVar,iBlockP,  i0, j1, k1]
    c110 = DA[iVar,iBlockP,  i1, j1, k0]
    c101 = DA[iVar,iBlockP,  i1, j0, k1]
    c111 = DA[iVar,iBlockP,  i1, j1, k1]

    c00 = c000*(1.-xd) + c100*xd
    c01 = c001*(1.-xd) + c101*xd
    c10 = c010*(1.-xd) + c110*xd
    c11 = c011*(1.-xd) + c111*xd

    c0 = c00*(1.-yd) + c10*yd
    c1 = c01*(1.-yd) + c11*yd

    c = c0*(1.-zd) + c1*zd
    return c

@njit
def interpolate_vectorized(points, iVar, DA, iTree_IA, node2block, pr):
    ret = np.empty(points.shape[0], dtype=np.float32)
    for i in range(points.shape[0]):
        ret[i] = interpolate_jit(points[i,:], iVar, DA, iTree_IA, node2block, pr)
    return ret

def interpolate(filetag, points, var='p', cache=None):
    """
    arguments:
        filetag:
            string with the name of the swmf output files,
            including the path, not includint the extension.
            e.g. '/home/user/data/3d__var_3_e20010101-010000-000'
        point(s):

        var (optional):
            string for the swmf variable name. Default: var='p'
    returns:
    """
    if cache is None:
        cache = read_all(filetag)
    elif filetag is not None:
        assert(cache['filetag'] == filetag)

    if len(points.shape) == 1:
        return    interpolate_jit(points, str2index[var], cache['DataArray'], cache['iTree_IA'], cache['node2block'], cache['pr'])
    return interpolate_vectorized(points, str2index[var], cache['DataArray'], cache['iTree_IA'], cache['node2block'], cache['pr'])


def B_dipole(x,y,z):
    DipoleStrength = 3.12e+4 #"dipole moment"(not really) in  nT * R_e**3  # https://en.wikipedia.org/wiki/Dipole_model_of_the_Earth%27s_magnetic_field
    M = np.array([0,0,DipoleStrength], dtype=np.float32)

    #ret = np.empty(X.shape, dtype=np.float32)
    #divr = 1./np.sqrt(X[:,0]**2 + X[:,1]**2 + X[:,2]**2)
    #ret[:,0] = ( 3.*(M[0]*X[:,0]+M[1]*X[:,1]+M[2]*X[:,2])*divr**5 )* X[:,0]  -  (divr**3)*M[0]
    #ret[:,1] = ( 3.*(M[0]*X[:,0]+M[1]*X[:,1]+M[2]*X[:,2])*divr**5 )* X[:,1]  -  (divr**3)*M[1]
    #ret[:,2] = ( 3.*(M[0]*X[:,0]+M[1]*X[:,1]+M[2]*X[:,2])*divr**5 )* X[:,2]  -  (divr**3)*M[2]
    #return ret

    divr = 1./np.sqrt(x**2 + y**2 + z**2)
    ret_x = ( 3.*(M[0]*x+M[1]*y+M[2]*z)*divr**5 )* x  -  (divr**3)*M[0]
    ret_y = ( 3.*(M[0]*x+M[1]*y+M[2]*z)*divr**5 )* y  -  (divr**3)*M[1]
    ret_z = ( 3.*(M[0]*x+M[1]*y+M[2]*z)*divr**5 )* z  -  (divr**3)*M[2]
    return ret_x, ret_y, ret_z


def swmf2vtk(filetag, use_ascii=False, cache=None):
    if cache is None:
        if filetag=='/tmp/3d__var_dipole':
            cache = read_all('/tmp/3d__var_2_e20190902-041000-000')
        else:
            cache = read_all(filetag)
    else:
        assert(cache['filetag'] == filetag)

    nBlock, nI, nJ, nK = cache['nBlock'], cache['nI'], cache['nJ'], cache['nK']

    x_blk = cache['DataArray'][_x,:,:,:,:]
    y_blk = cache['DataArray'][_y,:,:,:,:]
    z_blk = cache['DataArray'][_z,:,:,:,:]

    if filetag=='/tmp/3d__var_dipole':
        bx_blk, by_blk, bz_blk = B_dipole(x_blk.ravel(), y_blk.ravel(), z_blk.ravel())
        bx_blk = bx_blk.reshape((nBlock,nI,nJ,nK))
        by_blk = by_blk.reshape((nBlock,nI,nJ,nK))
        bz_blk = bz_blk.reshape((nBlock,nI,nJ,nK))
    else:
        bx_blk = cache['DataArray'][_bx,:,:,:,:]
        by_blk = cache['DataArray'][_by,:,:,:,:]
        bz_blk = cache['DataArray'][_bz,:,:,:,:]

    field = np.column_stack([bx_blk.ravel(),
                             by_blk.ravel(),
                             bz_blk.ravel()])

    npts = nBlock*nI*nJ*nK
    needed = nBlock*(nI+1)*(nJ+1)*(nK+1)

    all_vertices = np.empty( (needed, 3) )
    for iBlockP in range(nBlock):
        gridspacing = x_blk[iBlockP, 1,0,0] - x_blk[iBlockP, 0,0,0]

        xmin = x_blk[iBlockP, 0,0,0] - gridspacing/2.
        ymin = y_blk[iBlockP, 0,0,0] - gridspacing/2.
        zmin = z_blk[iBlockP, 0,0,0] - gridspacing/2.

        xmax = x_blk[iBlockP, nI-1,0   ,0   ] + gridspacing/2.
        ymax = y_blk[iBlockP, 0   ,nJ-1,0   ] + gridspacing/2.
        zmax = z_blk[iBlockP, 0   ,0   ,nK-1] + gridspacing/2.

        grid = np.mgrid[float(xmin):float(xmax+gridspacing):float(gridspacing), 
                        float(ymin):float(ymax+gridspacing):float(gridspacing),
                        float(zmin):float(zmax+gridspacing):float(gridspacing) ]
        grid = np.array(grid.reshape((3,(nI+1)*(nJ+1)*(nK+1))).transpose(), order='C')

        start = iBlockP*(nI+1)*(nJ+1)*(nK+1)
        end = (iBlockP+1)*(nI+1)*(nJ+1)*(nK+1)
        all_vertices[start:end,:] = grid

    unique_vertices, pointTo = np.unique(all_vertices,axis=0,return_inverse=True)
    assert(np.all( unique_vertices[pointTo, :] == all_vertices ))

    loc_in_all = np.arange(needed).reshape( (nBlock,(nI+1),(nJ+1),(nK+1)) )

    cells = []
    for iBlockP in range(nBlock):
        for i in range(nI):
            for j in range(nJ):
                for k in range(nK):
                    if True: # use vtk voxel for cells
                        celltype = 'VOXEL'
                        cells.append(
                             (pointTo[loc_in_all[iBlockP,i  ,j  ,k  ]] ,
                              pointTo[loc_in_all[iBlockP,i+1,j  ,k  ]] ,
                              pointTo[loc_in_all[iBlockP,i  ,j+1,k  ]] ,
                              pointTo[loc_in_all[iBlockP,i+1,j+1,k  ]] ,
                              pointTo[loc_in_all[iBlockP,i  ,j  ,k+1]] ,
                              pointTo[loc_in_all[iBlockP,i+1,j  ,k+1]] ,
                              pointTo[loc_in_all[iBlockP,i  ,j+1,k+1]] ,
                              pointTo[loc_in_all[iBlockP,i+1,j+1,k+1]] )
                            )
                    else: # use vtk hexahedron for cells
                        celltype = 'HEXAHEDRON'
                        cells.append(
                             (pointTo[loc_in_all[iBlockP,i  ,j  ,k  ]] ,
                              pointTo[loc_in_all[iBlockP,i+1,j  ,k  ]] ,
                              pointTo[loc_in_all[iBlockP,i+1,j+1,k  ]] ,
                              pointTo[loc_in_all[iBlockP,i  ,j+1,k  ]] ,
                              pointTo[loc_in_all[iBlockP,i  ,j  ,k+1]] ,
                              pointTo[loc_in_all[iBlockP,i+1,j  ,k+1]] ,
                              pointTo[loc_in_all[iBlockP,i+1,j+1,k+1]] ,
                              pointTo[loc_in_all[iBlockP,i  ,j+1,k+1]] )
                            )

    cells = np.array(cells, dtype=int)

    nVertices = unique_vertices.shape[0]
    nCells = npts

    if use_ascii:
        ftype='ASCII'
    else:
        ftype='BINARY'

    vtk_export(filetag+'.vtk', unique_vertices,
                    dataset = 'UNSTRUCTURED_GRID',
                    connectivity = {'CELLS' : {celltype : cells} },
                    cell_data = field,
                    texture = 'VECTORS',
                    cell_data_name = 'b',
                    ftype=ftype)


def swmf2vtk_separate(filetag, epsilon, use_ascii=False, cache=None):
    if cache is None:
        if filetag=='/tmp/3d__var_dipole':
            cache = read_all('/tmp/3d__var_2_e20190902-041000-000')
        else:
            cache = read_all(filetag)
    else:
        assert(cache['filetag'] == filetag)

    nBlock, nI, nJ, nK = cache['nBlock'], cache['nI'], cache['nJ'], cache['nK']

    cell_data_names = []
    textures = []
    cell_data = []

    x_blk = cache['DataArray'][_x,:,:,:,:]
    y_blk = cache['DataArray'][_y,:,:,:,:]
    z_blk = cache['DataArray'][_z,:,:,:,:]

    is_selected = np.empty(nBlock, dtype=bool)
    is_selected[:]=True
    is_selected = epsilon == x_blk[:, 1,0,0] - x_blk[:, 0,0,0]

    for vv in ['b','j','u','b1']:
        cell_data_names.append(vv)
        textures.append('VECTORS')
        cell_data.append(np.column_stack([cache['DataArray'][str2index[vv+'x'],is_selected,:,:,:].ravel(),
                                          cache['DataArray'][str2index[vv+'y'],is_selected,:,:,:].ravel(),
                                          cache['DataArray'][str2index[vv+'z'],is_selected,:,:,:].ravel()])
                        )
    for sv in ['rho','p']:
        cell_data_names.append(sv)
        textures.append('SCALARS')
        cell_data.append(cache['DataArray'][str2index[sv],is_selected,:,:,:].ravel())

    nSelected = np.count_nonzero(is_selected)

    all_vertices = []
    for iBlockP in range(nBlock):
        if not is_selected[iBlockP]: continue
        gridspacing = x_blk[iBlockP, 1,0,0] - x_blk[iBlockP, 0,0,0]

        xmin = x_blk[iBlockP, 0,0,0] - gridspacing/2.
        ymin = y_blk[iBlockP, 0,0,0] - gridspacing/2.
        zmin = z_blk[iBlockP, 0,0,0] - gridspacing/2.

        xmax = x_blk[iBlockP, nI-1,0   ,0   ] + gridspacing/2.
        ymax = y_blk[iBlockP, 0   ,nJ-1,0   ] + gridspacing/2.
        zmax = z_blk[iBlockP, 0   ,0   ,nK-1] + gridspacing/2.

        grid = np.mgrid[float(xmin):float(xmax+gridspacing):float(gridspacing), 
                        float(ymin):float(ymax+gridspacing):float(gridspacing),
                        float(zmin):float(zmax+gridspacing):float(gridspacing) ]
        grid = np.array(grid.reshape((3,(nI+1)*(nJ+1)*(nK+1))).transpose(), order='C')

        start = iBlockP*(nI+1)*(nJ+1)*(nK+1)
        end = (iBlockP+1)*(nI+1)*(nJ+1)*(nK+1)
        all_vertices.append(grid)
    all_vertices = np.vstack(all_vertices)

    unique_vertices, pointTo = np.unique(all_vertices,axis=0,return_inverse=True)
    assert(np.all( unique_vertices[pointTo, :] == all_vertices ))

    loc_in_block = np.arange((nI+1)*(nJ+1)*(nK+1)).reshape( ((nI+1),(nJ+1),(nK+1)) )
    cells = []

    startOfBlock = 0
    for iBlockP in range(nBlock):
        if not is_selected[iBlockP]: continue
        for i in range(nI):
            for j in range(nJ):
                for k in range(nK):
                    # use vtk voxel for cells
                    celltype = 'VOXEL'
                    cells.append(
                         (pointTo[startOfBlock+loc_in_block[i  ,j  ,k  ]] ,
                          pointTo[startOfBlock+loc_in_block[i+1,j  ,k  ]] ,
                          pointTo[startOfBlock+loc_in_block[i  ,j+1,k  ]] ,
                          pointTo[startOfBlock+loc_in_block[i+1,j+1,k  ]] ,
                          pointTo[startOfBlock+loc_in_block[i  ,j  ,k+1]] ,
                          pointTo[startOfBlock+loc_in_block[i+1,j  ,k+1]] ,
                          pointTo[startOfBlock+loc_in_block[i  ,j+1,k+1]] ,
                          pointTo[startOfBlock+loc_in_block[i+1,j+1,k+1]] )
                        )
        startOfBlock += (nI+1)*(nJ+1)*(nK+1)

    cells = np.array(cells, dtype=int)

    if use_ascii:
        ftype='ASCII'
    else:
        ftype='BINARY'

    outname = f'{filetag}_eps={epsilon}.vtk'
    outname = f'./eps={epsilon}.vtk'
    vtk_export(outname, unique_vertices,
                    dataset = 'UNSTRUCTURED_GRID',
                    connectivity = {'CELLS' : {celltype : cells} },
                    cell_data = cell_data,
                    texture = textures,
                    cell_data_name = cell_data_names,
                    ftype=ftype)


def model_api(filetag, cache=None):
    if cache is None:
        cache = read_all(filetag)
    else:
        assert(cache['filetag'] == filetag)

    #from magnetosphere import MagnetosphereState
    #MagnetosphereState(run_name=, time_Tstring=, interpolator=, data_array=, variables=, units=)


# python -c "from swmf_file_reader.read_swmf_files import swmf2vtk_separate; swmf2vtk_separate('',8.)"
if __name__ == '__main__':
    points = np.array([[-146.,  -14.,  -14.]])
    ftag = "/home/gary/temp/3d__var_3_e20031120-070000-000"
    #cac = read_all(ftag)
    #print(interpolate(ftag,point,var='p', cache=cac))
    #indx = find_index(ftag,point, cache=cac)
    #print(indx)
    #print(cac['DataArray'][(14,)+indx]) # (14,*indx) syntax only works python 3
    #read_all('/tmp/3d__var_dipole')
    #read_all('/tmp/3d__var_2_e20190902-041000-000')
    print(interpolate(ftag, points))
