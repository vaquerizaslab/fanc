from setuptools import setup, find_packages

__version__ = None
exec(open('kaic/version.py').read())

setup(
    name='kaic',
    version=__version__,
    description='Hi-C data analysis tools.',
    packages=find_packages(exclude=["test"]),
    install_requires=[
        'numpy',
        'scipy',
        'matplotlib',
        'pandas',
        'h5py',
        'pysam',
        'biopython',
        'pytest',
        'msgpack-python',
        'gridmap',
        'sklearn',
        'progressbar2',
        'pybedtools',
        'pyBigWig',
        'tables>=3.2.3'
    ],
    extras_require={
        'plotting': ["seaborn"]
    },
    scripts=['bin/kaic', 'bin/klot']
)
