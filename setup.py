from setuptools import setup, find_packages

setup(
    name='kaic',
    version='0.4.0',
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
        'tables>3.2.2'
    ],
    dependency_links=[
        "https://github.com/PyTables/PyTables/zipball/develop#egg=tables-3.3.0"
    ],
    extras_require={
        'plotting':  ["seaborn"]
    },
    scripts=['bin/kaic']
)
