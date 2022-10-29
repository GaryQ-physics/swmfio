import numpy as np
import scipy.io as sio

def read_batsrus(file):

    import os
    import swmfio


    (dirname, fname, fext) = swmfio.util.fileparts(file)
    assert fext == "" or fext == ".out" or fext == ".cdf"

    swmfio.logger.info("Creating class for file = " + file)
    if fext == '.cdf':
        from swmfio.batsrus_class import get_class_from_cdf
        cls = get_class_from_cdf(file)
        swmfio.logger.info("Created class for file = " + file)
        return cls 
    else:
        file = os.path.join(dirname, fname)
        from swmfio.batsrus_class import get_class_from_native
        cls = get_class_from_native(file)
        swmfio.logger.info("Created class for file = " + file)
        return cls


def read_tree(filetag):

    import swmfio

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
    nRoot_D = ff.read_ints(dtype=np.int32)  # The number of root nodes in all dimensions
    iTree_IA = ff.read_ints(dtype=np.int32).reshape((nInfo, nNode), order='F')
    ff.close()

    swmfio.logger.info(f"nDim = {nDim}")
    swmfio.logger.info(f"nInfo = {nInfo}")
    swmfio.logger.info(f"nNode = {nNode}")
    swmfio.logger.info(f"iRatio_D = {iRatio_D}")
    swmfio.logger.info(f"nRoot_D = {nRoot_D}")
    #swmfio.logger.info(f"iTree_IA = {iTree_IA}")

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

    import swmfio

    meta = {}
    nDim = 3
    ff = sio.FortranFile(filetag + ".out", 'r')

    header = ff.read_ints(dtype=np.uint8).tobytes().decode('UTF-8')
    nStep, Time, nDimOut, nParam, nVar = ff.read_ints(dtype=np.int32)

    n_D = ff.read_ints(dtype=np.int32)
    swmfio.logger.info(f"n_D = {n_D}")
    npts = n_D[0];
    assert(n_D[1] == 1 and n_D[2] == 1 and n_D.size==3)

    ScalarValues = ff.read_reals(dtype=np.float32)
    swmfio.logger.info(f"ScalarValues = {ScalarValues}")

    # nVar does not include x, y, and z.
    nVar = nVar + 3
    variables = ff.read_ints(dtype=np.uint8).tobytes().decode('UTF-8')
    variables = variables.strip().lower().split(' ')

    arrays = tuple(variables[:nVar]) # all other variables in the string are not in arrays in the file
    scalars = tuple(variables[nVar:])

    # +1 for variable 'measure' added later (volume)
    data = np.empty((npts, nVar+1), order='C', dtype=np.float32)
    data[:, -1] = np.nan # Should not be needed, but can be used to check.

    swmfio.logger.info(f"Reading {nVar} arrays")
    # Read xyz
    data[:,0:3] = ff.read_reals(dtype=np.float32).reshape(3, npts).T
    # Read grid variables
    for iVar in range(3, nVar):
        data[:, iVar] = ff.read_reals(dtype=np.float32)
    swmfio.logger.info(f"Read {nVar} arrays")

    try:
        ff.read_reals(dtype=np.float32)
        assert(False)
    except TypeError: # nessesary for python 2 compatibility, since it seems that there is no sio._fortran.FortranEOFError without it
        pass
    except sio._fortran.FortranEOFError:
        pass

    ff.close()

    meta['header'] = header.strip()
    meta['nStep'] = nStep
    meta['Time'] = Time
    meta['nDimOut'] = nDimOut
    meta['nParam'] = nParam
    meta['nVar'] = nVar
    meta['ScalarValues'] = ScalarValues
    meta['Arrays'] = arrays
    meta['Scalars'] = scalars

    units = meta['header'].split(' ')
    if units[0][0].isdigit():
        # In some files, the header contains a timestamp followed by units.
        # TODO: Find more general way to handle this. What is spec for header?
        meta['ArrayUnits'] = tuple(units[1:nVar+1])
    else:
        meta['ArrayUnits'] = tuple(units[:nVar])

    swmfio.logger.info("header: " + header.strip())
    swmfio.logger.info("arrays: {}".format(arrays))
    swmfio.logger.info("data.shape: {}".format(data.shape))

    return data, arrays, meta
