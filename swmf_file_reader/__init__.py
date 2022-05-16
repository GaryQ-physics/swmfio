def read_batsrus(file):
    from swmf_file_reader.batsrus_class import get_class_from_cdf
    from swmf_file_reader.batsrus_class import get_class_from_native

    if file.endswith('cdf'):
        return get_class_from_cdf(file)
    else:
        return get_class_from_native(file)


def read_rim(file):
    from swmf_file_reader.read_ie_files import read_iono_tec
    from swmf_file_reader.read_ie_files import read_iono_cdf

    if file.endswith('cdf'):
        return read_iono_cdf(file)
    else:
        return read_iono_tec(file)


def write_vtk(file, **kwargs):
    from swmf_file_reader import swmf2vtk
    swmf2vtk.write(file, **kwargs)