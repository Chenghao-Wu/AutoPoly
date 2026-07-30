from setuptools import setup, find_packages

VERSION = '2.0.0'
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
        package_data={'AutoPoly': ['extern/*']},
        include_package_data=True,
        packages=find_packages(),
        install_requires=['numpy >= 1.8.0',
                        'rdkit >= 2022.09.1'],
                        # add any additional packages that
        # needs to be installed along with your package. Eg: 'caer'

        keywords=['python', 'molecular dynamics'],
        extras_require={
            'dev': [
                'pytest>=7.0.0',
                'pytest-cov>=4.0.0',
                'pytest-mock>=3.10.0',
            ],
            'docs': [
                'mkdocs>=1.6.0',
                'mkdocs-material>=9.5.0',
                'mkdocstrings[python]>=0.26.0',
                'pymdown-extensions>=10.0',
            ],
        },
)