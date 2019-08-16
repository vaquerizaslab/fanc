import os
from setuptools import setup, find_packages, Command, Extension


__version__ = None
exec(open('kaic/version.py').read())


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
    name='kaic',
    version=__version__,
    description='Hi-C data analysis tools.',
    setup_requires=[
        'setuptools>=18.0',
        'cython'
    ],
    packages=find_packages(),
    package_data={'kaic': ['test/data/*/*']},
    install_requires=[
        'numpy>=1.16.0',
        'scipy',
        'pillow',
        'matplotlib',
        'pandas>=0.15.0',
        'pysam>=0.9.1',
        'biopython',
        'pytest',
        'msgpack-python',
        'msgpack-numpy',
        'gridmap',
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
        'genomic_regions>=0.0.4',
        'scikit-image>=0.15.0',
    ],
    extras_require={
        'cooler': ['cooler>=0.8.0', 'h5py'],
    },
    scripts=['bin/kaic', 'bin/klot'],
    cmdclass={
        'clean': CleanCommand
    },
    ext_modules=[
        Extension(
            'kaic.tools.sambam',
            sources=['kaic/tools/sambam.pyx', 'kaic/tools/natural_cmp.c'],
        ),
    ],
)
