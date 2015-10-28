from setuptools import setup, find_packages

setup(
    name='kaic',
    version='0.4.0',
    description='Hi-C data analysis tools.',
    packages=find_packages(exclude=["test"]),
    install_requires=[
        #'tables>=3.2.2',
        'pysam',
        'matplotlib',
        'pandas',
        'biopython',
        'numpy'
    ],
    extras_require={
        'plotting':  ["seaborn"]
    },
    scripts=['bin/bin_hic', 'bin/correct_hic', 'bin/filter_pairs', 'bin/filter_reads', 'bin/iterative_mapping',
             'bin/load_reads', 'bin/merge_hic', 'bin/mirny_to_kaic', 'bin/pairs_to_hic', 'bin/plot_hic_correlation',
             'bin/plot_hic_matrix', 'bin/plot_ligation_error', 'bin/reads_to_pairs']
)
