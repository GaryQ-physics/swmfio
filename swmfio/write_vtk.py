
def write_vtk(file_or_class, variables="all", epsilon=None, blocks=None, use_ascii=False):

    import numpy as np
    import swmfio

    variables_all = ['b', 'j', 'u', 'b1', 'rho', 'p', 'measure', 'block_id']
    if isinstance(variables, str):
        if variables == "all":
            variables = variables_all
            save_all = True
        else:
            variables = [variables]

    from swmfio import logger

    if isinstance(file_or_class, str):
        batsclass = swmfio.read_batsrus(file_or_class)
        fileout = file_or_class
    else:
        fileout = file_or_class.file
        batsclass = file_or_class

    logger.info("Creating VTK data structure.")

    DA = batsclass.DataArray
    vidx = batsclass.varidx

    nI = batsclass.nI
    nJ = batsclass.nJ
    nK = batsclass.nK
    nBlock = batsclass.block2node.size
    nNode = batsclass.node2block.size

    logger.debug(" (nI, nJ, nK) = ({0:d}, {1:d}, {2:d})".format(nI, nJ, nK))
    logger.debug(" nBlock = {0:d}".format(nBlock))
    logger.debug(" nNode  = {0:d}".format(nNode))

    nVar = len(batsclass.varidx)

    assert(DA.shape == (nVar, nI, nJ, nK, nBlock))
    assert(np.isfortran(DA))

    x_blk = DA[vidx['x'],:,:,:,:]
    y_blk = DA[vidx['y'],:,:,:,:]
    z_blk = DA[vidx['z'],:,:,:,:]
    
    is_selected = np.full(nBlock, True, dtype=bool)

    if epsilon is not None:
        is_selected = epsilon == x_blk[1,0,0, :] - x_blk[0,0,0, :]

    if blocks is not None:
        blocks = np.array(blocks)
        is_selected[:] = False
        is_selected[blocks] = True

    cell_data = []
    for vv in ['b','j','u','b1']:
        if not vv in variables: continue
        cell_data.append(
            {
                "name" : vv,
                "texture" : "VECTORS",
                "array" : np.column_stack([DA[vidx[vv+'x'],:,:,:,is_selected].ravel(),
                                           DA[vidx[vv+'y'],:,:,:,is_selected].ravel(),
                                           DA[vidx[vv+'z'],:,:,:,is_selected].ravel()])
            })
    for sv in ['rho','p', 'measure']:
        if not sv in variables: continue
        cell_data.append(
            {
                "name" : sv,
                "texture" : "SCALARS",
                "array" : DA[vidx[sv],:,:,:,is_selected].ravel()
            })

    nSelected = np.count_nonzero(is_selected)

    logger.debug(" nSelected = {0:d} (of {1:d})".format(nSelected, nBlock))

    block_id = np.full(nSelected*(nI)*(nJ)*(nK), -1, dtype=np.int32)
    all_vertices = np.full((nSelected*(nI+1)*(nJ+1)*(nK+1), 3), np.nan, dtype=np.float32)
    startOfBlock = 0    # Counter of points in block
    cellIndexStart = 0  # Counter of cells in block
    logger.info(" Creating block grids.")
    for iBlockP in range(nBlock):

        logger.debug(f"  Creating grid for block #{iBlockP+1}/{nBlock+1}")

        if blocks is not None and iBlockP > blocks[-1]:
            logger.debug("  iBlockP > blocks[-1]. Done.")
            break

        if not is_selected[iBlockP]:
            logger.debug(f"  Block #{iBlockP+1} not selected. Omitting.")
            continue

        block_id[cellIndexStart:cellIndexStart+nI*nJ*nK] = iBlockP
        cellIndexStart = cellIndexStart + nI*nJ*nK

        gridspacing = x_blk[1,0,0, iBlockP] - x_blk[0,0,0, iBlockP]

        xmin = x_blk[0,0,0, iBlockP] - gridspacing/2.
        ymin = y_blk[0,0,0, iBlockP] - gridspacing/2.
        zmin = z_blk[0,0,0, iBlockP] - gridspacing/2.

        xmax = x_blk[nI-1,0   ,0   , iBlockP] + gridspacing/2.
        ymax = y_blk[0   ,nJ-1,0   , iBlockP] + gridspacing/2.
        zmax = z_blk[0   ,0   ,nK-1, iBlockP] + gridspacing/2.

        logger.debug("    (x, y, z) min = ({0:.1f}, {1:.1f}, {2:.1f})".format(xmin, ymin, zmin))
        logger.debug("    (x, y, z) max = ({0:.1f}, {1:.1f}, {2:.1f})".format(xmax, ymax, zmax))

        grid = np.mgrid[float(xmin):float(xmax+gridspacing):float(gridspacing),
                        float(ymin):float(ymax+gridspacing):float(gridspacing),
                        float(zmin):float(zmax+gridspacing):float(gridspacing) ]
        grid = np.array(grid.reshape((3,(nI+1)*(nJ+1)*(nK+1))).transpose(), order='C')

        all_vertices[startOfBlock:startOfBlock+(nI+1)*(nJ+1)*(nK+1), :] = grid
        startOfBlock += (nI+1)*(nJ+1)*(nK+1)

    logger.info(" Created block grids.")


    unique_vertices, pointTo = np.unique(all_vertices, axis=0, return_inverse=True)
    assert(np.all( unique_vertices[pointTo, :] == all_vertices ))

    # Nodes of blocks
    loc_in_block = np.arange((nI+1)*(nJ+1)*(nK+1)).reshape( ((nI+1),(nJ+1),(nK+1)) )
    cells = []
    startOfBlock = 0

    celltype = 'VOXEL' # or HEXAHEDRON
    logger.info(f' Creating {celltype}s.')
    for iBlockP in range(nBlock):

        logger.debug(f"  Creating cells for block #{iBlockP+1}/{nBlock+1}")

        if blocks is not None and iBlockP > blocks[-1]:
            logger.debug("  iBlockP > blocks[-1]. Done.")
            break

        if not is_selected[iBlockP]:
            logger.debug(f"  Block #{iBlockP+1} not selected. Omitting.")
            continue

        # TODO: These loops can be vectorized.
        for i in range(nI):
            for j in range(nJ):
                for k in range(nK):

                    if celltype == 'VOXEL':
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
                    if celltype == 'HEXAHEDRON':
                        cells.append(
                             (pointTo[startOfBlock+loc_in_block[i  ,j  ,k  ]] ,
                              pointTo[startOfBlock+loc_in_block[i+1,j  ,k  ]] ,
                              pointTo[startOfBlock+loc_in_block[i+1,j+1,k  ]] ,
                              pointTo[startOfBlock+loc_in_block[i  ,j+1,k  ]] ,
                              pointTo[startOfBlock+loc_in_block[i  ,j  ,k+1]] ,
                              pointTo[startOfBlock+loc_in_block[i+1,j  ,k+1]] ,
                              pointTo[startOfBlock+loc_in_block[i+1,j+1,k+1]] ,
                              pointTo[startOfBlock+loc_in_block[i  ,j+1,k+1]] )
                            )

        startOfBlock += (nI+1)*(nJ+1)*(nK+1)


    cells = np.array(cells, dtype=int)

    logger.info(f' Created {celltype}s.')    

    logger.info("Created VTK data structure.")

    if 'block_id' in variables:
        cell_data.append(
            {
                "name" : "block_id",
                "texture" : "SCALARS",
                "array" : block_id
            })

    if use_ascii:
        ftype='ASCII'
    else:
        ftype='BINARY'

    extra = ""
    if epsilon is not None:
        extra = f"_epsilon={epsilon}"
    if variables is not None and save_all == False:
        extra = "_vars=" + ",".join(variables)
    if blocks is not None:
        if len(blocks) < 5:
            blockstrs = []
            for block in blocks:
                blockstrs.append(str(block))
            extra = "_blocks=" + ",".join(blockstrs)
        else:
            extra = f"_Nblocks={len(blocks)}"
    fileout = fileout + extra + ".vtk"

    debug = False
    if logger.getEffectiveLevel() > 20:
        debug = True

    logger.info("Writing " + fileout)

    from swmfio.vtk_export import vtk_export
    vtk_export(fileout, unique_vertices,
                    dataset='UNSTRUCTURED_GRID',
                    connectivity={'CELLS': {celltype: cells} },
                    cell_data=cell_data,
                    ftype=ftype,
                    debug=debug)

    logger.info("Wrote " + fileout)

    return fileout
