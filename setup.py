from setuptools import setup, find_packages

setup(
    name='kaic',
    version='0.4.0',
    description='Hi-C data analysis tools.',
    packages=find_packages(exclude=["test"]),
    install_requires=[
        'tables>=3.2.2',
        'pysam',
        'matplotlib',
        'pandas',
        'biopython',
        'numpy'
    ],
    extras_require={
        'plotting':  ["seaborn"]
    }
)
