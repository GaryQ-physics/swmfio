import numpy as np
from numba import njit
import re

def grep_dash_o(RE, lines):
    ret = ''
    for line in lines:
        findall = re.findall(RE, line)
        if findall != []:
            ret = ret + '\n'.join(findall) + '\n'
    return ret

@njit
def unravel_index(index, shape, order='C'):
    if order=='F':
        pass
    else:
        assert(order=='C')
        shape = shape[::-1]

    multiindex = np.empty(len(shape), dtype=np.int64)

    linearind = index
    for dim in range(len(shape)):
        quotient, remainder = divmod(linearind, shape[dim])
        multiindex[dim] = remainder
        linearind = quotient

    if order == 'F':
        return multiindex
    else:
        return multiindex[::-1]
