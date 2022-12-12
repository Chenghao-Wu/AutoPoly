from setuptools import setup, find_packages

VERSION = '0.0.1' 
DESCRIPTION = 'AutoPoly'
LONG_DESCRIPTION = 'Build your LAMMPS data file'

# Setting up
setup(
       # the name must match the folder name 'verysimplemodule'
        name="AutoPoly", 
        version=VERSION,
        author="Zhenghao Wu",
        author_email="zhenghao.wu95@gmail.com",
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        packages=find_packages(),
        install_requires=['numpy >= 1.8.0',
                        'rdkit >= 2022.09.1',
                        'pathlib >= 1.0.1'], 
                        # add any additional packages that 
        # needs to be installed along with your package. Eg: 'caer'
        
        keywords=['python', 'molecular dynamics'],
)