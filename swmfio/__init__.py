from swmfio.read_rim import read_rim
from swmfio.read_batsrus import read_batsrus
from swmfio.write_vtk import write_vtk

def logger():
    import sys
    import logging
    logger = logging.getLogger('swmfio')
    handler = logging.StreamHandler(stream=sys.stdout)
    formatter = logging.Formatter('%(asctime)s.%(msecs)03d:%(filename)s:%(funcName)s(): %(message)s', datefmt='%H:%M:%S')
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    for handler in logger.handlers[:]:
        logger.removeHandler(handler)
    logger.addHandler(handler)
    logger.setLevel(logging.INFO)

    return logger

logger = logger()

