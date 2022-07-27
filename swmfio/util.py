import numpy as np
from numba import njit
import re

def _dlfile(file, tmpdir=None):

    import os
    import swmfio

    if tmpdir is None:
        import tempfile
        import platform
        system = platform.system()
        if system in ["Darwin", "Linux"]:
            tmpdir = "/tmp"
        else:
            tmpdir = tempfile.gettempdir()

    (dirname, fname, fext) = fileparts(file)
    filename = fname + fext
    subdir = dirname.replace("http://","").replace("https://","")
    subdir = os.path.join(tmpdir, subdir)
    if not os.path.exists(subdir):
        swmfio.logger.info("Creating " + subdir)
        os.makedirs(subdir)
    tmppath = os.path.join(subdir, filename)

    if os.path.exists(tmppath):
        swmfio.logger.info("Found " + tmppath)
    else:
        partfile = tmppath + ".part"
        swmfio.logger.info("Downloading\n  " + file + "\n  to\n  " + partfile)
        try:
            import urllib.request
            #urllib.request.urlretrieve(dirname + "/" + filename, partfile)
            req = urllib.request.Request(dirname + "/" + filename)
            response = urllib.request.urlopen(req)
            #print(response.info().get('Content-Encoding'))
            length = response.info().get('Content-Length')
            length_downloaded = 0
            CHUNK = 16 * 1024
            import sys
            with open(partfile, 'wb') as f:
              while True:
                if length is not None:
                    length_downloaded = length_downloaded + CHUNK
                    sys.stdout.write(7*"\b" + "{0:.2f}{1:s}".format(100*length_downloaded/int(length), "%"))
                    sys.stdout.flush()
                chunk = response.read(CHUNK)
                if not chunk:
                  sys.stdout.write("\n")  
                  break
                f.write(chunk)

        except Exception as e:
            swmfio.logger.info("Download error")
            # Won't remove .part file if application crashes.
            # Would need to register an atexit.
            if os.path.exists(partfile):
                swmfio.logger.info("Removing " + partfile)
                os.remove(partfile)
            raise 

        swmfio.logger.info("Renaming\n  " + partfile + "\n  to  \n  " + tmppath)
        os.rename(partfile, tmppath)

    return tmppath

def dlfile(file, tmpdir=None):

    import os
    import swmfio

    (dirname, fname, fext) = swmfio.util.fileparts(file)

    if fext == "" or fext == ".out":
        swmfio.logger.info("Remote files are raw simulation output.")
        for ext in ['.tree', '.info', '.out']:
            tmpfile = _dlfile(file + ext, tmpdir=tmpdir)
        (dirname, fname, fext) = swmfio.util.fileparts(tmpfile)
        tmpfile = os.path.join(dirname, fname)         
    else:
        swmfio.logger.info("Remote file is CCMC CDF.")
        tmpfile = _dlfile(file, tmpdir=tmpdir)        
        tmpfile = os.path.join(dirname, tmpfile)

    return tmpfile

def fileparts(file):
    import os
    (dirname, fname) = os.path.split(file)
    (fname, fext) = os.path.splitext(fname)

    if file.startswith("http"):
        fname_split = file.split("/")
        fname, fext = os.path.splitext(fname_split[-1])
        dirname = "/".join(fname_split[0:-1])
    else:
        file_parts = fname

    return dirname, fname, fext

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
