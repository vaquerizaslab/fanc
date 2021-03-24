import os
from setuptools import setup, find_packages, Command, Extension


__version__ = None
exec(open('fanc/version.py').read())


class CleanCommand(Command):
    """
    Custom clean command to tidy up the project root.
    """
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        os.system('rm -vrf ./build ./dist ./*.pyc ./*.tgz ./*.egg-info ./htmlcov')


setup(
    name='fanc',
    version=__version__,
    description='Framework for the ANalysis of C-data.',
    setup_requires=[
        'setuptools>=18.0',
        'cython'
    ],
    packages=find_packages(),
    package_data={'fanc': ['test/data/*/*']},
    include_package_data=True,
    install_requires=[
        'numpy>=1.16.0',
        'scipy',
        'pillow',
        'matplotlib>=3.1.0',
        'pandas>=0.15.0',
        'pysam>=0.9.1',
        'biopython',
        'pytest',
        'msgpack>=1.0.0',
        'msgpack-numpy>=0.4.6.1',
        'scikit-learn',
        'progressbar2',
        'pybedtools',
        'pyBigWig',
        'PyYAML>=5.1',
        'tables>=3.5.1',
        'seaborn',
        'future',
        'gridmap>=0.14.0',
        'intervaltree',
        'genomic_regions>=0.0.10',
        'scikit-image>=0.15.0',
        'cooler>=0.8.0',
        'h5py'
    ],
    scripts=['bin/fanc', 'bin/fancplot'],
    cmdclass={
        'clean': CleanCommand
    },
    ext_modules=[
        Extension(
            'fanc.tools.sambam',
            sources=['fanc/tools/sambam.pyx', 'fanc/tools/natural_cmp.c'],
        ),
    ],
)
