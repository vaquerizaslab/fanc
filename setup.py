import os
from setuptools import setup, find_packages, Command

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
    # packages=find_packages(exclude=["test"]),
    install_requires=[
        'numpy>=1.8.0',
        'scipy',
        'matplotlib',
        'pandas>=0.15.0',
        'pysam',
        'biopython',
        'pytest',
        'msgpack-python',
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
    }
)
