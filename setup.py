from setuptools import setup, find_packages

install_requires = ["numpy","scipy",
                    # numba, magnetovis@''
                    ]

setup(
    name='swmf_file_reader',
    version='0.0.1',
    packages=find_packages(),
    description='reading swmf files using spacepy and scipy.io.FortranFiles',
    install_requires=install_requires
     )
