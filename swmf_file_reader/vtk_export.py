import os
import numpy as np

def cleanUponError(writefile):
    def clean_writefile(*args, **kwargs):
        if os.path.exists(args[0]) and ('debug' in kwargs.keys()) and kwargs['debug']:
            print(f'warning: overwiting existing {args[0]}')

        try:
            writefile(*args, **kwargs)
        except:
            os.remove(args[0])
            raise Exception (f'removing bad file {args[0]}')

    return clean_writefile

@cleanUponError
def vtk_export(out_filename, points,
                    dataset = 'UNSTRUCTURED_GRID',
                    connectivity = None,
                    point_data = None,
                    cell_data = None,
                    title='Title',
                    ftype='BINARY',
                    debug=True):
    """
    Parameters
    ----------
    out_filename : string
        the full filename (including path) that the vtk will be saved as.
    points : np.ndarray
        Nx3 numpy array where points[n,i] is the i'th cartesian component of the nth point
    dataset : string, optional
        the type of dataset for vtk (e.g. 'STRUCTURED_GRID' or 'POLYDATA').
        The default is 'UNSTRUCTURED_GRID'.
    connectivity : dictionary, optional
        details all connectivity types involved with this dataset. The default is {} (empty dictionary).
    point_data : list of dictionaries, optional
        The default is None.
        keys : "array", "name", "texture"
    cell_data: list of dictionaries, optional
        The default is None.
        keys : "array", "name", "texture"
    title : string, optional
        The default is 'Title'.
    ftype : string, optional
        The default is 'BINARY'.

    Note : some of the above arguments can instead be keyword strings that implement common behavior.
           e.g.: dataset='POLYDATA', connectivity='LINES'
 
    Returns
    -------
    None.

    point_data = [{"array": ..., "name":"rho", "texture":"SCALARS"}, {"array": ..., "name":"b", "texture":"VECTORS"}]
    """
    if connectivity is None:
        connectivity = {}

    if isinstance(point_data, dict):
        point_data = [point_data,]

    if isinstance(cell_data, dict):
        cell_data = [cell_data,]

    if dataset == 'POLYDATA' and connectivity == 'LINES':
        connectivity = {'LINES' : np.array([points.shape[0]])}

    if dataset == 'STRUCTURED_GRID' and type(connectivity) != dict:
        connectivity = {'DIMENSIONS' : tuple(connectivity)}

    f = open(out_filename,'wb')
    if debug:
        print("Writing " + out_filename)
    f.write(b'# vtk DataFile Version 3.0\n')
    f.write(b'%s\n'%(bytearray(title, 'utf-8')))
    f.write(b'%s\n'%(bytearray(ftype, 'utf-8')))
    f.write(b'DATASET %s\n'%(bytearray(dataset, 'utf-8')))

    num_points = points.shape[0]
    if dataset == 'STRUCTURED_GRID':
        Nx,Ny,Nz = connectivity['DIMENSIONS']
        f.write(b'DIMENSIONS %d %d %d\n'%(Nx, Ny, Nz))
        assert(Nx*Ny*Nz == num_points)

    f.write(b'POINTS %d float\n'%(num_points))

    if ftype=='BINARY':
        points = np.array(points, dtype='>f')
        f.write(points.tobytes())
    elif ftype=='ASCII':
        np.savetxt(f, points)

    f.write(b'\n')

    for key in connectivity.keys():
        if key == 'CELLS':
            celltypes = {
                            'VERTEX'         : 1 ,
                            'POLY_VERTEX'    : 2 ,
                            'LINE'           : 3 ,
                            'POLY_LINE'      : 4 ,
                            'TRIANGLE'       : 5 ,
                            'TRIANGLE_STRIP' : 6 ,
                            'POLYGON'        : 7 ,
                            'PIXEL'          : 8 ,
                            'QUAD'           : 9 ,
                            'TETRA'          : 10,
                            'VOXEL'          : 11,
                            'HEXAHEDRON'     : 12,
                            'WEDGE'          : 13,
                            'PYRAMID'        : 14
                        }
            ctypes = list(connectivity['CELLS'].keys())
            for i in range(len(ctypes)):
                ctype = ctypes[i]
                cint = celltypes[ctype]
                cnum = connectivity['CELLS'][ctype].shape[0]
                cnpts = connectivity['CELLS'][ctype].shape[1]
                if i==0:
                    carr = np.column_stack([ cnpts*np.ones(cnum, dtype=int),
                                             connectivity['CELLS'][ctype]
                                           ])
                    ct_arr = cint*np.ones(cnum, dtype=int)
                else:
                    carr = np.vstack([ carr,
                                       np.column_stack([ cnpts*np.ones(cnum, dtype=int),
                                                         connectivity['CELLS'][ctype]
                                                       ])
                                     ])
                    ct_arr.append(cint*np.ones(cnum, dtype=int))

            num_cells = ct_arr.size
            f.write(b'CELLS %d %d\n'%(num_cells, carr.size))
            if ftype=='BINARY':
                carr = np.array(carr, dtype=">i4")
                f.write(carr.tobytes())
            else:
                np.savetxt(f,carr,fmt='%d')
            f.write(b'CELL_TYPES %d\n'%(ct_arr.size))
            if ftype=='BINARY':
                ct_arr = np.array(ct_arr, dtype=">i4")
                f.write(ct_arr.tobytes())
            else:
                np.savetxt(f,ct_arr,fmt='%d')

        if key == 'LINES':
            num_curves = len(connectivity['LINES'])
            f.write(b'LINES %d %d\n'%(num_curves, num_points + num_curves))
            lines = np.zeros(num_points+num_curves, dtype=int)
            i = 0
            j = 0
            for header in connectivity['LINES']:
                lines[i] = header
                lines[i+1 : i+header+1] = np.arange(j, j + header)
                i = i + header + 1
                j = j + header

            if ftype=='BINARY':
                lines = np.array(lines, dtype=">i4")
                f.write(lines.tobytes())
            elif ftype=='ASCII':
                np.savetxt(f, lines, fmt='%d')

    def writedata(name, array, texture): #todo integer values??
        if ' ' in name:
            raise ValueError ('names of point and cell data cannot contain spaces')

        if texture == 'SCALARS':
            f.write(b'SCALARS %s float 1\n'%(bytearray(name, 'utf-8'))) # number with float???
            f.write(b'LOOKUP_TABLE default\n')
        if texture == 'VECTORS':
            f.write(b'VECTORS %s float\n'%(bytearray(name, 'utf-8')))
        if texture == 'TEXTURE_COORDINATES':
            f.write(b'TEXTURE_COORDINATES %s 2 float\n'%(bytearray(name, 'utf-8'))) # http://www.earthmodels.org/data-and-tools/topography/paraview-topography-by-texture-mapping
        if ftype=='BINARY':
            array = np.array(array, dtype='>f')
            f.write(array.tobytes())
        elif ftype=='ASCII':
            np.savetxt(f, array)

    if point_data is not None:
        f.write(b'POINT_DATA %d\n'%(num_points))
        for data in point_data:
            writedata(data['name'], data['array'], data['texture'])

    if cell_data is not None:
        f.write(b'CELL_DATA %d\n'%(num_cells))
        for data in cell_data:
            writedata(data['name'], data['array'], data['texture'])

    f.close()
    if debug:
        print("Wrote " + out_filename)
        print("Open in ParaView on command line using magnetovis.sh --data=" + out_filename)
