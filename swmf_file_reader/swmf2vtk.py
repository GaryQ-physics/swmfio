import numpy as np
import swmf_file_reader.batsrus_class as batscls
from swmf_file_reader.vtk_export import vtk_export

def write(filetag, debug=False, epsilon=None, use_ascii=False):

    if isinstance(filetag, str): # TODO: Check extenstion
        if debug:
            print("Reading {}.*".format(filetag))
        bats_slice = batscls.get_class_from_native(filetag)
        if debug:
            print("Read {}.*".format(filetag))
        outname = filetag
    else:
        bats_slice = filetag
        outname = '/tmp/write_BATSRUS_unstructured_grid_vtk'

    if debug:
        print("Creating VTK data structure")

    DA = bats_slice.DataArray
    vidx = bats_slice.varidx

    nI = bats_slice.nI
    nJ = bats_slice.nJ
    nK = bats_slice.nK
    nBlock = bats_slice.block2node.size
    nVar = len(bats_slice.varidx)

    assert(DA.shape == (nVar, nI, nJ, nK, nBlock))
    assert(np.isfortran(DA))

    x_blk = DA[vidx['x'],:,:,:,:]
    y_blk = DA[vidx['y'],:,:,:,:]
    z_blk = DA[vidx['z'],:,:,:,:]

    is_selected = np.empty(nBlock, dtype=bool)
    is_selected[:] = True
    if epsilon is not None:
        is_selected = epsilon == x_blk[1,0,0, :] - x_blk[0,0,0, :]

    cell_data = []
    for vv in ['b','j','u','b1']:
        cell_data.append({
            "name" : vv,
            "texture" : "VECTORS",
            "array" : np.column_stack([DA[vidx[vv+'x'],:,:,:,is_selected].ravel(),
                                       DA[vidx[vv+'y'],:,:,:,is_selected].ravel(),
                                       DA[vidx[vv+'z'],:,:,:,is_selected].ravel()])
                        })
    for sv in ['rho','p', 'measure']:
        cell_data.append({
            "name" : sv,
            "texture" : "SCALARS",
            "array" : DA[vidx[sv],:,:,:,is_selected].ravel()
                        })

    nSelected = np.count_nonzero(is_selected)
    all_vertices = np.empty((nSelected*(nI+1)*(nJ+1)*(nK+1), 3), dtype=np.float32)
    all_vertices[:,:] = np.nan
    startOfBlock = 0
    for iBlockP in range(nBlock):
        if not is_selected[iBlockP]: continue
        gridspacing = x_blk[1,0,0, iBlockP] - x_blk[0,0,0, iBlockP]

        xmin = x_blk[0,0,0, iBlockP] - gridspacing/2.
        ymin = y_blk[0,0,0, iBlockP] - gridspacing/2.
        zmin = z_blk[0,0,0, iBlockP] - gridspacing/2.

        xmax = x_blk[nI-1,0   ,0   , iBlockP] + gridspacing/2.
        ymax = y_blk[0   ,nJ-1,0   , iBlockP] + gridspacing/2.
        zmax = z_blk[0   ,0   ,nK-1, iBlockP] + gridspacing/2.

        grid = np.mgrid[float(xmin):float(xmax+gridspacing):float(gridspacing), 
                        float(ymin):float(ymax+gridspacing):float(gridspacing),
                        float(zmin):float(zmax+gridspacing):float(gridspacing) ]
        grid = np.array(grid.reshape((3,(nI+1)*(nJ+1)*(nK+1))).transpose(), order='C')

        all_vertices[startOfBlock:startOfBlock+(nI+1)*(nJ+1)*(nK+1), :] = grid
        startOfBlock += (nI+1)*(nJ+1)*(nK+1)

    #print(is_selected)
    #print(np.all(is_selected))
    #print(all_vertices)
    #print(np.count_nonzero(np.isnan(all_vertices)))

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
                    if True: # use vtk voxel for cells
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
                    else: # use vtk hexahedron for cells
                        celltype = 'HEXAHEDRON'
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

    if debug:
        print("Created VTK data structure")

    if use_ascii:
        ftype='ASCII'
    else:
        ftype='BINARY'

    if epsilon is not None:
        outname = f'{outname}_epsilon={epsilon}.vtk'
    else:
        outname = f'{outname}.vtk'

    vtk_export(outname, unique_vertices,
                    dataset='UNSTRUCTURED_GRID',
                    connectivity={'CELLS': {celltype: cells} },
                    cell_data=cell_data,
                    ftype=ftype,
                    debug=debug)


if __name__ == '__main__':
    write('/tmp/3d__var_2_e20190902-041000-000', epsilon=None)
