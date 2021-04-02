import numpy as np
from swmf_constants import Used_,Status_,Level_,Parent_,Child0_,Child1_,Coord1_,CoordLast_,ROOTNODE_
from vtk_export_copy import vtk_export

def F2P(fortran_index):
    return fortran_index - 1

def P2F(python_index):
    return python_index + 1


def read_info_file(filetag):
    dINFO={}
    with open(filetag+'.info','r') as f:
        for line in f.readlines():
            if line == '\n' : continue
            if line[0] == '#': continue
            splt = line.split()
            if len(splt) == 2:
                dINFO[splt[1]] = splt[0]
    return dINFO


def read_tree_file(filetag):
    # first read info file
    dINFO=read_info_file(filetag)
    nDim = int(dINFO['nDim'])
    nI = int(dINFO['BlockSize1'])
    nJ = int(dINFO['BlockSize2'])
    nK = int(dINFO['BlockSize3'])

    ## Loading AMR tree
    # directly read bytes
    #f = open(filetag+".tree", 'rb')
    #filebytes = np.fromfile(f, dtype=np.int32)
    #f.close()

    try:
        # use scipy FortranFile
        from scipy.io import FortranFile
        ff = FortranFile(filetag+".tree", 'r')
        if True:
            nDim, nInfo, nNode = ff.read_reals(dtype=np.int32)
            iRatio_D = ff.read_reals(dtype=np.int32) # Array of refinement ratios
            nRoot_D = ff.read_reals(dtype=np.int32) # The number of root nodes in all dimension
            iTree_IA = ff.read_reals(dtype=np.int32).reshape((nInfo,nNode), order='F')
        else:
            nDim, nInfo, nNode = ff.read_ints(dtype='i4')
            iRatio_D = ff.read_ints(dtype='i4') # Array of refinement ratios
            nRoot_D = ff.read_ints(dtype='i4') # The number of root nodes in all dimension
            iTree_IA = ff.read_ints(dtype='i4').reshape((nInfo,nNode), order='fortran')
    except:
        raise RuntimeWarning ("scipy.io.FortranFile didnt work")
        # use fortranfile
        from fortranfile import FortranFile
        ff = FortranFile(filetag+".tree") # read or write ???
        nDim, nInfo, nNode = ff.readInts()
        iRatio_D = ff.readInts() # Array of refinement ratios
        nRoot_D = ff.readInts() # The number of root nodes in all dimension
        iTree_IA = ff.readInts().reshape((nInfo,nNode), order='fortran')

    ########################### check_thing_work #######################
    # Maximum number of ghost cells set by Config.pl script.
    # Valid values are 0,1,2,3,4,5
    nG = 2
    # Refinement ratios in the 3 dimensions. Either 1 or 2.
    # The values are set by the Config.pl script.
    iRatio, jRatio, kRatio = min(2, nI), min(2, nJ), min(2, nK)
    # Number of dimensions in which grid adaptation is done
    nDimAmr = iRatio + jRatio + kRatio - 3
    assert(nDimAmr == nDim)
    assert(nDim == 3)
    # Number of children per node
    nChild = 2**nDimAmr
    assert(np.isfortran(iTree_IA))
    assert(np.all(iRatio_D == np.array([iRatio, jRatio, kRatio])))
    ####################################################################

    dTREE = {'iTree_IA': iTree_IA, 'nRoot_D': nRoot_D, 'iRatio_D': iRatio_D,
             'nDims': 3, 'nI': nI, 'nJ': nJ, 'nK': nK }
    dTREE['xGlobalMin'] = float(dINFO['Coord1Min'])
    dTREE['yGlobalMin'] = float(dINFO['Coord2Min'])
    dTREE['zGlobalMin'] = float(dINFO['Coord3Min'])
    dTREE['xGlobalMax'] = float(dINFO['Coord1Max'])
    dTREE['yGlobalMax'] = float(dINFO['Coord2Max'])
    dTREE['zGlobalMax'] = float(dINFO['Coord3Max'])

    return dTREE


def read_out_file(filetag):
    pass


# from SWMF/GM/BATSRUS/srcBATL/BATL_tree.f90 line 951 with substitutions
def get_tree_position(iNode, dTREE, returnall=False):
    '''
    Calculate normalized position of the edges of node iNode.
    Zero is at the minimum boundary of the grid, one is at the max boundary
    '''
    iTree_IA = dTREE['iTree_IA']
    nRoot_D = dTREE['nRoot_D']
    iRatio_D = dTREE['iRatio_D']

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

    if returnall:
        return PositionMin_D, PositionMax_D, MaxIndex_D, BlockCoord_D
    else:
        return PositionMin_D, PositionMax_D # what was returned in original


def get_physical_dimensions(iNode, dTREE, returnCenters=False):
    x_start = dTREE['xGlobalMin']
    y_start = dTREE['yGlobalMin']
    z_start = dTREE['zGlobalMin']
    x_range = dTREE['xGlobalMax'] - dTREE['xGlobalMin']
    y_range = dTREE['yGlobalMax'] - dTREE['yGlobalMin']
    z_range = dTREE['zGlobalMax'] - dTREE['zGlobalMin']

    iLevel = dTREE['iTree_IA'][F2P(Level_), F2P(iNode)]
    assert(dTREE['nI'] == dTREE['nJ'] == dTREE['nK'])
    assert(x_range == y_range == z_range)
    gridspacing = (x_range/dTREE['nI'])*0.5**iLevel

    PositionMin_D, PositionMax_D = get_tree_position(iNode, dTREE)
    xmin = x_range*(PositionMin_D[0]) + x_start
    ymin = y_range*(PositionMin_D[1]) + y_start
    zmin = z_range*(PositionMin_D[2]) + z_start
    xmax = x_range*(PositionMax_D[0]) + x_start
    ymax = y_range*(PositionMax_D[1]) + y_start
    zmax = z_range*(PositionMax_D[2]) + z_start

    xlims = (xmin+gridspacing/2., xmax-gridspacing/2.)
    ylims = (ymin+gridspacing/2., ymax-gridspacing/2.)
    zlims = (zmin+gridspacing/2., zmax-gridspacing/2.)

    xminmax = (xmin, xmax)
    yminmax = (ymin, ymax)
    zminmax = (zmin, zmax)

    if returnCenters:
        return xlims, ylims, zlims, gridspacing
    else:
        return xminmax, yminmax, zminmax, gridspacing


# supposed to reproduce SWMF/GM/BATSRUS/srcBATL/BATL_tree.f90 line 975, but differently
def find_tree_node(point, dTREE):
    iTree_IA = dTREE['iTree_IA']

    xin = dTREE['xGlobalMin'] <= point[0] <= dTREE['xGlobalMax']
    yin = dTREE['yGlobalMin'] <= point[1] <= dTREE['yGlobalMax']
    zin = dTREE['zGlobalMin'] <= point[2] <= dTREE['zGlobalMax']

    if not (xin and yin and zin): 
        raise RuntimeError ('point out of simulation volume')

    iNode = ROOTNODE_
    while(True):
        if Used_ == iTree_IA[F2P(Status_), F2P(iNode)]:
            break

        for j in range(8):
            child = iTree_IA[F2P(Child1_+j), F2P(iNode)]
            xminmax, yminmax, zminmax, gridspacing = get_physical_dimensions(child, dTREE, returnCenters=False)

            xin = xminmax[0] <= point[0] <= xminmax[1]
            yin = yminmax[0] <= point[1] <= yminmax[1]
            zin = zminmax[0] <= point[2] <= zminmax[1]

            if xin and yin and zin: 
                iNode = child
                break

    return iNode


def read_all(filetag):
    def Reorder(arr1, arr2):
        ''' arr1 and arr2 are (N,3) arrays of for N  3d-points 
            if arr1 and arr2 contain the same points, but in different
            order, returns array of indices 'ind' such that
            arr1[ind,:] == arr2
        '''
        arr1 = np.array(arr1); arr2 = np.array(arr2)
        assert(len(arr1.shape) == len(arr2.shape) == 2)
        if arr1.shape != arr2.shape: return False

        sort1 = np.lexsort((arr1[:,0],arr1[:,1],arr1[:,2]))  
        sort2 = np.lexsort((arr2[:,0],arr2[:,1],arr2[:,2]))
        #undo1 = np.argsort(sort1)
        undo2 = np.argsort(sort2)

        if not np.all(arr2 == (arr1[sort1[undo2],:])):
            raise RuntimeError ('arrays inputed arent reordering of each other')
        return sort1[undo2]

    if True:
        import spacepy.pybats.bats as bats
        import read_swmf_files as rswmf
        data = bats.Bats2d(filetag + ".out")
        header = "R R R Mp/cc km/s km/s km/s J/m3 nT nT nT nT nT nT nPa uA/m2 uA/m2 uA/m2 --"
        assert(header == data.meta['header'].strip())
    else:
        data = read_out_file(filetag) # TODO : eventually
    dTREE = read_tree_file(filetag)

    # in what follows:
    #  the P in iNodeP and iBlockP stands for python like indexing (as oposed to fortran)
    #  
    #  iNodeP indexes all nodes of the tree, from 0 to nNode-1,
    #  and thus the "iNode" to be used in the other functions is simply iNodeP+1, or P2F(iNodeP)
    # 
    #  iBlockP indexes all the blocks used, from 0 to nBlock-1.
    #  There is one for each node with a status of used. 
    #  Note, nBlock*nI*nJ*nK = total number of batsrus cells (npts)
    points = np.column_stack([data['x'],data['y'],data['z']])
    npts = points.shape[0]
    nI, nJ, nK = dTREE['nI'], dTREE['nJ'], dTREE['nK']
    nBlock = npts//(nI*nJ*nK) if npts%(nI*nJ*nK)==0 else -1
    nNode = dTREE['iTree_IA'].shape[1]

    iBlockP = 0
    block2node = -np.ones((nBlock,), dtype=int)
    node2block = -np.ones((nNode,), dtype=int)
    reconstructed_points = np.nan*np.empty((npts,3), dtype=np.float32)
    for iNodeP in range(nNode):
        if dTREE['iTree_IA'][F2P(Status_), iNodeP] == Used_:
            block2node[iBlockP] = iNodeP
            node2block[iNodeP] = iBlockP
            xlims, ylims, zlims, gridspacing = rswmf.get_physical_dimensions(P2F(iNodeP), dTREE, returnCenters=True)
            assert(xlims[1]-xlims[0] == (nI-1)*gridspacing)
            assert(ylims[1]-ylims[0] == (nJ-1)*gridspacing)
            assert(zlims[1]-zlims[0] == (nK-1)*gridspacing)

            grid = np.mgrid[xlims[0]:xlims[1]+gridspacing:gridspacing, 
                            ylims[0]:ylims[1]+gridspacing:gridspacing,
                            zlims[0]:zlims[1]+gridspacing:gridspacing ]
            grid = np.array(grid.reshape((3,nI*nJ*nK)).transpose(), order='C')

            #### equivalent too ####
            #grid = np.empty((nI,nJ,nK,3))
            #for i in range(nI):
            #    for j in range(nJ):
            #        for k in range(nK):
            #            grid[i,j,k, 0] = xlims[0]+gridspacing*i
            #            grid[i,j,k, 1] = ylims[0]+gridspacing*j
            #            grid[i,j,k, 2] = zlims[0]+gridspacing*k
            #grid = grid.reshape((nI*nJ*nK,3))
            ############

            start = iBlockP*nI*nJ*nK
            end = (iBlockP+1)*nI*nJ*nK
            reconstructed_points[start:end,:] = grid

            iBlockP += 1

    ind = Reorder(points, reconstructed_points)
    #assert(np.all( points[ind,:] == reconstructed_points ))
    return [data, dTREE, ind, block2node, node2block]

def get_block_data(filetag):
    data, dTREE, ind, block2node, node2block = read_all(filetag)
    nBlock, nI, nJ, nK = block2node.size, dTREE['nI'], dTREE['nJ'], dTREE['nK']
    block_data = {}
    for key in 'x y z jx jy jz bx by bz b1x b1y b1z'.split(' '):
        block_data[key] = data[key][ind].reshape((nBlock, nI, nJ, nK))
    return [block_data, nBlock, nI, nJ, nK]

def find_index(point):
    pass


def interpolate(filetag, point, var='p', debug=False):
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
    #TODO: vectorize find_tree_node for points (N,3) (possibly with jitFORTRAN).
    #      then get rid of get_physical_dimensions() call in this function
    #      gridpacing and minmax can be found direcly from block_data.
    #      then maybe this function can be vectorized as well for points (N,3)
    # this function should maybe go in different file also

    data, dTREE, ind, block2node, node2block = read_all(filetag)
    nBlock, nI, nJ, nK = block2node.size, dTREE['nI'], dTREE['nJ'], dTREE['nK']

    def getvar(_var, iNode, i,j,k):
        return data[_var][ind].reshape(nBlock,nI,nJ,nK)[node2block[F2P(iNode)],i,j,k]

    iNode = find_tree_node(point, dTREE)

    # get the gridspacing in x,y,z
    gridspacingX = getvar('x', iNode, 1,0,0) - getvar('x', iNode, 0,0,0)
    gridspacingY = getvar('y', iNode, 0,1,0) - getvar('y', iNode, 0,0,0)
    gridspacingZ = getvar('z', iNode, 0,0,1) - getvar('z', iNode, 0,0,0)

    # i0 is s.t. the highest index s.t. the x coordinate of the 
    #  corresponding cell block_data[iNode,i0,:,:]  is still less than point[0]
    i0 = int(np.floor( (point[0] - getvar('x', iNode, 0, 0, 0))/gridspacingX ))
    j0 = int(np.floor( (point[1] - getvar('y', iNode, 0, 0, 0))/gridspacingY ))
    k0 = int(np.floor( (point[2] - getvar('z', iNode, 0, 0, 0))/gridspacingZ ))
    if debug: print(getvar('x', iNode, 0, 0, 0),getvar('y', iNode, 0, 0, 0),getvar('z', iNode, 0, 0, 0))
    if debug: print(i0,j0,k0)
    #i1 = i0+1 is the lowest index s.t. the x coordinate of the 
    #  corresponding cell block_data[iNode,i1,:,:]  is still greater than point[0]

    # together, i0 and i1 form the upper and lower bounds for a linear interpolation in x
    # likewise for j0,j1,y  and k0,k1,z

    #TODO implement better interpolation at ends of block.
    # This method effectively makes it nearest neighbor at the ends
    if i0 == -1:
        i0 = 0
        i1 = 0
        if debug: print('edge case')
    elif i0 == dTREE['nI']-1:
        i1 = dTREE['nI']-1
        if debug: print('edge case')
    else:
        i1 = i0 + 1

    if j0 == -1:
        j0 = 0
        j1 = 0
        if debug: print('edge case')
    elif j0 == dTREE['nJ']-1:
        j1 = dTREE['nJ']-1
        if debug: print('edge case')
    else:
        j1 = j0 + 1

    if k0 == -1:
        k0 = 0
        k1 = 0
        if debug: print('edge case')
    elif k0 == dTREE['nK']-1:
        k1 = dTREE['nK']-1
        if debug: print('edge case')
    else:
        k1 = k0 + 1

    # all together i0,i1,j0,ect... form a cube of side length "gridpacing" to do trililear interpolation within
    # define xd as the distance along x of point within that cube, in units of "gridspacing"
    xd = (point[0] - getvar('x', iNode, i0, 0 , 0 ) )/gridspacingX
    yd = (point[1] - getvar('y', iNode, 0 , j0, 0 ) )/gridspacingY
    zd = (point[2] - getvar('z', iNode, 0 , 0 , k0) )/gridspacingZ

    if debug: print(xd,yd,zd)
    if debug: print(getvar('x', iNode,  i0, j0, k0))
    if debug: print(getvar('y', iNode,  i0, j0, k0))
    if debug: print(getvar('z', iNode,  i0, j0, k0))
    if debug: print('hellothere')
    if debug: print((iNode, node2block[F2P(iNode)], i0, j0, k0))


    #https://en.wikipedia.org/wiki/Trilinear_interpolation
    c000 = getvar(var, iNode,  i0, j0, k0)
    c001 = getvar(var, iNode,  i0, j0, k1)
    c010 = getvar(var, iNode,  i0, j1, k0)
    c100 = getvar(var, iNode,  i1, j0, k0)
    c011 = getvar(var, iNode,  i0, j1, k1)
    c110 = getvar(var, iNode,  i1, j1, k0)
    c101 = getvar(var, iNode,  i1, j0, k1)
    c111 = getvar(var, iNode,  i1, j1, k1)
    if debug: print(c000)
    if debug: print(c001,c010,c100,c100,c011,c110,c101,c111)

    c00 = c000*(1.-xd) + c100*xd
    c01 = c001*(1.-xd) + c101*xd
    c10 = c010*(1.-xd) + c110*xd
    c11 = c011*(1.-xd) + c111*xd

    c0 = c00*(1.-yd) + c10*yd
    c1 = c01*(1.-yd) + c11*yd

    c = c0*(1.-zd) + c1*zd
    if debug: print(c)
    return c


def swmf2vtk_old(filetag, use_ascii=False): #use vtk hexahedron for cells
    block_data, nBlock, nI, nJ, nK = get_block_data(filetag)

    npts = nBlock*nI*nJ*nK
    needed = nBlock*(nI+1)*(nJ+1)*(nK+1)

    all_vertices = np.empty( (needed, 3) )
    for iBlockP in range(nBlock):
        gridspacing = block_data['x'][iBlockP, 1,0,0] - block_data['x'][iBlockP, 0,0,0]

        xmin = block_data['x'][iBlockP, 0,0,0] - gridspacing/2.
        ymin = block_data['y'][iBlockP, 0,0,0] - gridspacing/2.
        zmin = block_data['z'][iBlockP, 0,0,0] - gridspacing/2.

        xmax = block_data['x'][iBlockP, nI-1,0   ,0   ] + gridspacing/2.
        ymax = block_data['y'][iBlockP, 0   ,nJ-1,0   ] + gridspacing/2.
        zmax = block_data['z'][iBlockP, 0   ,0   ,nK-1] + gridspacing/2.

        grid = np.mgrid[xmin:xmax+gridspacing:gridspacing, 
                        ymin:ymax+gridspacing:gridspacing,
                        zmin:zmax+gridspacing:gridspacing ]
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

    field = np.column_stack([block_data['bx'].ravel(),
                             block_data['by'].ravel(),
                             block_data['bz'].ravel()])

    if use_ascii:
        ftype='ASCII'
    else:
        ftype='BINARY'
    vtk_export(filetag+'-hexahedron.vtk', unique_vertices,
                    dataset = 'UNSTRUCTURED_GRID',
                    connectivity = {'CELLS' : {'HEXAHEDRON': cells} },
                    cell_data = field,
                    texture = 'VECTORS',
                    cell_data_name = 'b',
                    ftype=ftype)



def swmf2vtk(filetag, use_ascii=False): #use vtk voxel for cells
    block_data, nBlock, nI, nJ, nK = get_block_data(filetag)

    npts = nBlock*nI*nJ*nK
    needed = nBlock*(nI+1)*(nJ+1)*(nK+1)

    all_vertices = np.empty( (needed, 3) )
    for iBlockP in range(nBlock):
        gridspacing = block_data['x'][iBlockP, 1,0,0] - block_data['x'][iBlockP, 0,0,0]

        xmin = block_data['x'][iBlockP, 0,0,0] - gridspacing/2.
        ymin = block_data['y'][iBlockP, 0,0,0] - gridspacing/2.
        zmin = block_data['z'][iBlockP, 0,0,0] - gridspacing/2.

        xmax = block_data['x'][iBlockP, nI-1,0   ,0   ] + gridspacing/2.
        ymax = block_data['y'][iBlockP, 0   ,nJ-1,0   ] + gridspacing/2.
        zmax = block_data['z'][iBlockP, 0   ,0   ,nK-1] + gridspacing/2.

        grid = np.mgrid[xmin:xmax+gridspacing:gridspacing, 
                        ymin:ymax+gridspacing:gridspacing,
                        zmin:zmax+gridspacing:gridspacing ]
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
    cells = np.array(cells, dtype=int)

    nVertices = unique_vertices.shape[0]
    nCells = npts

    field = np.column_stack([block_data['bx'].ravel(),
                             block_data['by'].ravel(),
                             block_data['bz'].ravel()])

    if use_ascii:
        ftype='ASCII'
    else:
        ftype='BINARY'
    vtk_export(filetag+'.vtk', unique_vertices,
                    dataset = 'UNSTRUCTURED_GRID',
                    connectivity = {'CELLS' : {'VOXEL': cells} },
                    cell_data = field,
                    texture = 'VECTORS',
                    cell_data_name = 'b',
                    ftype=ftype)



if __name__ == '__main__':
    print(interpolate("/home/gary/temp/3d__var_3_e20031120-070000-000",point,var='p', debug=True))

