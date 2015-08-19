from setuptools import setup, find_packages

setup(
    name='kaic',
    version='0.3.0',
    description=('Vacuum Hi-C data analysis tools.'),
    packages = find_packages(exclude=["test"])
)