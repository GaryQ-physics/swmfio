from setuptools import setup, find_packages

install_requires = ["numpy","scipy","numba","cdflib>=0.4.4"]

setup(
    name='swmf_file_reader',
    version='0.0.4',
    author='Gary Quaresima, Bob Weigel',
    author_email='garyquaresima@gmail.com,rweigel@gmu.edu',
    packages=find_packages(),
    description='Fast read and interpolation of SWMF/BATSRUS native .out and CCMC .cdf files. Also reads SWMF/RIM .tec and .cdf files.',
    install_requires=install_requires
)
