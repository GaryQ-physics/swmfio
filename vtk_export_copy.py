import numpy as np


def vtk_export(out_filename, points,
                    dataset = 'UNSTRUCTURED_GRID',
                    connectivity = {},
                    point_data = None,
                    cell_data = None,
                    texture = None,
                    point_data_name = 'point_data',
                    cell_data_name = 'cell_data',
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
    point_data : TYPE, optional
        DESCRIPTION. The default is None.
    texture : string, optional
        
    point_data_name : TYPE, optional
        DESCRIPTION. The default is 'point data'.
    title : TYPE, optional
        DESCRIPTION. The default is 'Title'.
    ftype : TYPE, optional
        DESCRIPTION. The default is 'BINARY'.

    Returns
    -------
    None.
    """

    #if grid == 'STRUCTURED_GRID':
    #    Nx,Ny,Nz = connectivity

    if dataset == 'POLYDATA' and connectivity == 'LINES':
        connectivity = {'LINES' : np.zeros(points.shape[0])}

    if dataset == 'STRUCTURED_GRID' and type(connectivity) != dict:
        connectivity = {'DIMENSIONS' : tuple(connectivity)}

    #assert(out_filename[0] == '/') #!!!!!!!!!!!!!!!!!
    f = open(out_filename,'w')
    if debug:
        print("Writing " + out_filename)
    f.write('# vtk DataFile Version 3.0\n')
    f.write(title + '\n')
    f.write(ftype + '\n')
    f.write('DATASET ' + dataset + '\n')

    num_points = points.shape[0]
    if dataset == 'STRUCTURED_GRID':
        Nx,Ny,Nz = connectivity['DIMENSIONS']
        f.write('DIMENSIONS ' + str(Nx) + ' ' + str(Ny) + ' ' + str(Nz) + '\n' )
        assert(Nx*Ny*Nz == num_points)
        #f.write('POINTS '+str(Nx*Ny*Nz)+' float\n')
    #if dataset == 'POLYDATA' or dataset == 'UNSTRUCTURED_GRID':
        #f.write('\nPOINTS '+str(num_points)+' float\n')

    f.write('POINTS ' + str(num_points) + ' float\n')

    if ftype=='BINARY':
        points = np.array(points, dtype='>f')
        f.write(points.tobytes())
    elif ftype=='ASCII':
        np.savetxt(f, points)

    f.write('\n')

    for key in connectivity.keys():
        #if key == 'TRIANGLE_STRIPS':
        if key == 'CELLS':
            celltypes = {#TODO fix whitespace
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
            for i in range(len(connectivity['CELLS'].keys())):
                ctype = connectivity['CELLS'].keys()[i]
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

            f.write('CELLS %d %d\n'%(ct_arr.size,carr.size))
            if ftype=='BINARY':
                carr = np.array(carr, dtype=">i4")
                f.write(carr.tobytes())
            else:
                np.savetxt(f,carr,fmt='%d')
            f.write('CELL_TYPES %d\n'%(ct_arr.size))
            if ftype=='BINARY':
                ct_arr = np.array(ct_arr, dtype=">i4")
                f.write(ct_arr.tobytes())
            else:
                np.savetxt(f,ct_arr,fmt='%d')

        if key == 'LINES':
            num_curves = len(connectivity['LINES'])
            f.write("LINES {} {} \n".format(num_curves, num_points + num_curves))
            lines = np.zeros(num_points+num_curves)
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

        if key == 'HYPERBOLOID TRIANGLE':
            assert(ftype == 'ASCII')
            extra_circle = 1
            closing_points = 2
            head_count = 1
            extra_point = 1
            points_on_circle = 0
            points_on_circle = connectivity[key]
            
            num_circles = int(num_points/points_on_circle)
            
            triangle_line = ''
            lines_of_triangles = num_circles - extra_circle
            head_triangle = 2 * points_on_circle + closing_points
            points_on_triangle_strips = (2 
                                        * points_on_circle 
                                        + closing_points 
                                        + head_count) * lines_of_triangles
            f.write('\n')
            f.write('TRIANGLE_STRIPS {} {}\n'.format(lines_of_triangles, 
                                                     points_on_triangle_strips))
            for i in range(num_points - points_on_circle):
                triangle_line = triangle_line + ' {} {}'.format(i, i+points_on_circle)
                
                if (i+1) % int(points_on_circle) == 0:
                    if i == points_on_circle: 
                        f_repeat_point = 0
                        l_repeat_point = f_repeat_point + points_on_circle
                    else: 
                        f_repeat_point = i - points_on_circle + extra_point
                        l_repeat_point = f_repeat_point + points_on_circle
                    line = '{} {} {} {}\n'.format(head_triangle, 
                                               triangle_line,
                                               f_repeat_point,
                                               l_repeat_point)
                    f.write(line)
                    triangle_line = ''
            f.write('\n')

    if point_data is not None:
        if isinstance(point_data_name, list):

            f.write('POINT_DATA ' + str(num_points) + '\n')
            for i in range(len(point_data_name)):
                point_data_name_ELEMENT = point_data_name[i]
                point_data_ELEMENT = point_data[i]
                texture_ELEMENT = texture[i]

                if texture_ELEMENT == 'SCALARS':
                    f.write('SCALARS ' + point_data_name_ELEMENT + ' float 1\n') # number with float???
                    f.write('LOOKUP_TABLE default\n')
                if texture_ELEMENT == 'VECTORS':
                    f.write('VECTORS ' + point_data_name_ELEMENT + ' float\n')
                if texture_ELEMENT == 'TEXTURE_COORDINATES':
                    f.write('TEXTURE_COORDINATES ' + point_data_name_ELEMENT + ' 2 float\n') # http://www.earthmodels.org/data-and-tools/topography/paraview-topography-by-texture-mapping

                if ftype=='BINARY':
                    point_data_ELEMENT = np.array(point_data_ELEMENT, dtype='>f')
                    f.write(point_data_ELEMENT.tobytes())
                elif ftype=='ASCII':
                    np.savetxt(f, point_data_ELEMENT)

        else:
            f.write('POINT_DATA ' + str(num_points) + '\n')

            if texture == 'SCALARS':
                f.write('SCALARS ' + point_data_name + ' float 1\n') # number with float???
                f.write('LOOKUP_TABLE default\n')
            if texture == 'VECTORS':
                f.write('VECTORS ' + point_data_name + ' float\n')
            if texture == 'TEXTURE_COORDINATES':
                f.write('TEXTURE_COORDINATES ' + point_data_name + ' 2 float\n') # http://www.earthmodels.org/data-and-tools/topography/paraview-topography-by-texture-mapping

            if ftype=='BINARY':
                point_data = np.array(point_data, dtype='>f')
                f.write(point_data.tobytes())
            elif ftype=='ASCII':
                np.savetxt(f, point_data)

    if cell_data is not None:
        f.write('CELL_DATA %d\n'%(cell_data.shape[0]))
        if texture == 'VECTORS':
            f.write('VECTORS '+cell_data_name+' float\n')
        else:
            raise ValueError ('TODO implement')
        if ftype=='BINARY':
            cell_data = np.array(cell_data, dtype='>f')
            f.write(cell_data.tobytes())
        else:
            np.savetxt(f, cell_data, fmt = '%.3f')

    f.close()
    if debug:
        print("Wrote " + out_filename)
        print("Open in ParaView on command line using magnetovis.sh --data=" + out_filename)

