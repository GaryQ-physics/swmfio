import numpy as np
import scipy.io as sio

def read_batsrus(file):

    import os
    import swmfio

    swmfio.logger.info("Called with file = " + file)

    (dirname, fname, fext) = swmfio.util.fileparts(file)
    assert fext == "" or fext == ".out" or fext == ".cdf"

    if fext == '.cdf':
        from swmfio.batsrus_class import get_class_from_cdf
        return get_class_from_cdf(file)
    else:
        file = os.path.join(dirname, fname)
        from swmfio.batsrus_class import get_class_from_native
        return get_class_from_native(file)


def read_tree(filetag):
    # first read info file
    info = {'filetag' : filetag}
    with open(filetag+'.info','r') as f:
        for line in f.readlines():
            if line == '\n' : continue
            if line[0] == '#': continue
            splt = line.split()
            if len(splt) == 2:
                info[splt[1]] = splt[0]

    ## load tree file
    ff = sio.FortranFile(filetag+".tree", 'r')
    nDim, nInfo, nNode = ff.read_ints(dtype=np.int32)
    iRatio_D = ff.read_ints(dtype=np.int32) # Array of refinement ratios
    nRoot_D = ff.read_ints(dtype=np.int32) # The number of root nodes in all dimension
    iTree_IA = ff.read_ints(dtype=np.int32).reshape((nInfo,nNode), order='F')
    ff.close()

    assert(nDim == int(info['nDim']))

    return iTree_IA, iRatio_D, nRoot_D, info

    # ########################### check_things_work #######################
    # # Maximum number of ghost cells set by Config.pl script.
    # # Valid values are 0,1,2,3,4,5
    # nG = 2
    # # Refinement ratios in the 3 dimensions. Either 1 or 2.
    # # The values are set by the Config.pl script.
    # iRatio, jRatio, kRatio = min(2, pr.nI), min(2, pr.nJ), min(2, pr.nK)
    # # Number of dimensions in which grid adaptation is done
    # nDimAmr = iRatio + jRatio + kRatio - 3
    # assert(nDimAmr == nDim)
    # assert(nDim == 3)
    # # Number of children per node
    # nChild = 2**nDimAmr
    # assert(np.isfortran(iTree_IA))
    # assert(np.all(iRatio_D == np.array([iRatio, jRatio, kRatio])))
    # ####################################################################


def read_data(filetag):
    nDim = 3
    ff = sio.FortranFile(filetag+".out", 'r')
    header = ff.read_ints(dtype=np.uint8).tobytes().decode('UTF-8')
    nStep, Time, nDimOut, nParam, nVar_tmp = ff.read_ints(dtype=np.int32)
    n_D = ff.read_ints(dtype=np.int32)
    ScalarParams = ff.read_reals(dtype=np.float32)
    npts = n_D[0]; assert(n_D[1]==1 and n_D[2]==1 and n_D.size==3)

    # x,y,and z  are stored as one larger array in file. All others are stored seperately
    # as a results xyz count only once towards nVar.
    nVar = nVar_tmp + 2; del nVar_tmp
    variables = ff.read_ints(dtype=np.uint8).tobytes().decode('UTF-8')
    variables = variables.strip().split(' ')
    variables = tuple(variables[:nVar]) # all other variables in the string arent in arrays in the file

    data_arr = np.empty((npts, nVar+1), order='C', dtype=np.float32)
    data_arr[:, -1] = np.nan

    xyz = ff.read_reals(dtype=np.float32).reshape(3, npts) 
    data_arr[:, 0] = xyz[0,:]
    data_arr[:, 1] = xyz[1,:]
    data_arr[:, 2] = xyz[2,:]
    for iVar in range(3,nVar):
        data_arr[:, iVar] = ff.read_reals(dtype=np.float32)

    status = ff.read_reals(dtype=np.float32) #? whatever this means
    # check to make sure ff is at the end of the file
    try:
        ff.read_reals(dtype=np.float32)
        assert(False)
    except TypeError: # nessesary for python 2 compatibility, since it seems that there is no sio._fortran.FortranEOFError without it
        pass
    except sio._fortran.FortranEOFError:
        pass

    ff.close()
    expectedheader = "R R R Mp/cc km/s km/s km/s J/m3 nT nT nT nT nT nT nPa uA/m2 uA/m2 uA/m2 --"
    assert(expectedheader == header.strip())
    return data_arr, variables

