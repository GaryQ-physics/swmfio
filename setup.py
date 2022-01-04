from setuptools import setup, find_packages

install_requires = ["numpy","scipy","numba"]
# install_requires.append("cdflib")

setup(
    name='swmf_file_reader',
    version='0.0.3',
    author='Gary Quaresima, Bob Weigel',
    author_email='garyquaresima@gmail.com,rweigel@gmu.edu',
    packages=find_packages(),
    description='Fast read and interpolate of native SWMF files',
    install_requires=install_requires
)
