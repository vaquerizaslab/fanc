from setuptools import setup, find_packages

setup(
    name='kaic',
    version='0.1.1',
    description=('Vacuum Hi-C data analysis tools.'),
    packages = find_packages(exclude="test"),
    install_requires=[
        'hiclib',
        'mirnylib',
        'numpy',
        'matplotlib'
    ]
)