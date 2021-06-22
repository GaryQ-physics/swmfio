from setuptools import setup, find_packages

install_requires = ["numpy","scipy"]

setup(
    name='swmf_file_reader',
    version='0.0.0.8',
    packages=find_packages(),
    description='reading swmf files using spacepy and scipy.io.FortranFiles',
    install_requires=install_requires
     )
