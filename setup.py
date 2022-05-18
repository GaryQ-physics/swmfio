from setuptools import setup, find_packages

install_requires = ["numpy","scipy","numba","cdflib>=0.4.4"]

setup(
    name='swmfio',
    version='0.9.0',
    author='Gary Quaresima, Bob Weigel',
    author_email='garyquaresima@gmail.com,rweigel@gmu.edu',
    packages=find_packages(),
    description='Fast read and interpolation of SWMF/BATSRUS native .out and CCMC .cdf files. Also reads SWMF/RIM .tec and .cdf files.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    install_requires=install_requires
)
