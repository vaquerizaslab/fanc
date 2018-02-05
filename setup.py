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
        'numpy>=1.8.0',
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
        'PyYAML',
        'tables>=3.2.3',
        'seaborn',
        'future',
        'gridmap'
    ],
    scripts=['bin/kaic', 'bin/klot'],
    cmdclass={
        'clean': CleanCommand
    },
    ext_modules=[
        Extension(
            'sambam',
            sources=['kaic/tools/sambam.pyx'],
        ),
    ],
)
