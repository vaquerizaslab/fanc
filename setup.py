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
        'seaborn'
    ],
    scripts=['bin/kaic', 'bin/klot']
)
