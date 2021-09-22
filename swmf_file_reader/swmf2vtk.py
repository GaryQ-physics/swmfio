import numpy as np
import swmf_file_reader.batsrus_class as batscls
from magnetovis.vtk_export import vtk_export
from hxform import hxform as hx


def B_dipole(X):
    DipoleStrength = 3.12e+4 #"dipole moment"(not really) in  nT * R_e**3  # https://en.wikipedia.org/wiki/Dipole_model_of_the_Earth%27s_magnetic_field
    M = np.array([0,0,DipoleStrength], dtype=np.float32)

    #divr = 1./np.sqrt(x**2 + y**2 + z**2)
    #ret_x = ( 3.*(M[0]*x+M[1]*y+M[2]*z)*divr**5 )* x  -  (divr**3)*M[0]
    #ret_y = ( 3.*(M[0]*x+M[1]*y+M[2]*z)*divr**5 )* y  -  (divr**3)*M[1]
    #ret_z = ( 3.*(M[0]*x+M[1]*y+M[2]*z)*divr**5 )* z  -  (divr**3)*M[2]
    #return ret_x, ret_y, ret_z

    ret = np.empty(X.shape, dtype=np.float32)
    divr = 1./np.sqrt(X[:,0]**2 + X[:,1]**2 + X[:,2]**2)
    ret[:,0] = ( 3.*(M[0]*X[:,0]+M[1]*X[:,1]+M[2]*X[:,2])*divr**5 )* X[:,0]  -  (divr**3)*M[0]
    ret[:,1] = ( 3.*(M[0]*X[:,0]+M[1]*X[:,1]+M[2]*X[:,2])*divr**5 )* X[:,1]  -  (divr**3)*M[1]
    ret[:,2] = ( 3.*(M[0]*X[:,0]+M[1]*X[:,1]+M[2]*X[:,2])*divr**5 )* X[:,2]  -  (divr**3)*M[2]
    return ret


def get_dipole_field(xyz):
    # Xyz_D and returned b_D in SMG (SM) coordinates
    b = np.empty(xyz.shape, dtype=np.float32)
    r = np.sqrt(np.sum(xyz, axis=1))
    DipoleStrength = 3.12e+4 #"dipole moment"(not really) in  nT * R_e**3  # https://en.wikipedia.org/wiki/Dipole_model_of_the_Earth%27s_magnetic_field
    Term1      = DipoleStrength*xyz[:,2]*3/r**2
    b[:, 0:2] = Term1[:,None]*xyz[:, 0:2]/r[:,None]**3
    b[:, 2]    = (Term1*xyz[:,2]-DipoleStrength)/r**3
    return b


def write_BATSRUS_unstructured_grid_vtk(filetag, epsilon=None, use_ascii=False):
    if isinstance(filetag, str):# todo check extenstion
        bats_slice = batscls.get_class_from_native(filetag)
        outname = filetag
    else:
        bats_slice = filetag
        outname = '/tmp/write_BATSRUS_unstructured_grid_vtk'

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

    ####################################################################
    time = (2019, 9, 2, 4, 10, 0)
    #xyz_GSM = np.column_stack([x_blk.ravel(), y_blk.ravel(), z_blk.ravel()])
    xyz_GSM = np.column_stack([DA[vidx['x'],:,:,:,is_selected].ravel(),
                               DA[vidx['y'],:,:,:,is_selected].ravel(),
                               DA[vidx['z'],:,:,:,is_selected].ravel()])
    xyz_SMG = hx.transform(xyz_GSM, time, 'GSM', 'SM')
    dip_SMG = B_dipole(xyz_SMG)
    dip_GSM = hx.transform(dip_SMG, time, 'SM', 'GSM')

    #print(np.count_nonzero(np.isnan(xyz_GSM)))
    #print(np.count_nonzero(np.isnan(xyz_SMG)))
    #print(np.count_nonzero(np.isnan(dip_SMG)))
    #print(np.count_nonzero(np.isnan(dip_GSM)))

    cell_data.append({
        "name" : 'dipole',
        "texture" : "VECTORS",
        "array" : dip_GSM })
    del time
    ####################################################################

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

    if use_ascii:
        ftype='ASCII'
    else:
        ftype='BINARY'

    if epsilon is not None:
        outname = f'{outname}_epsilon={epsilon}.vtk'
    else:
        outname = f'{outname}.vtk'

    vtk_export(outname, unique_vertices,
                    dataset = 'UNSTRUCTURED_GRID',
                    connectivity = {'CELLS' : {celltype : cells} },
                    cell_data = cell_data,
                    ftype=ftype)


if __name__ == '__main__':
    write_BATSRUS_unstructured_grid_vtk('/tmp/3d__var_2_e20190902-041000-000', epsilon=None)
