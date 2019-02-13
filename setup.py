from setuptools import setup, find_packages

setup(
    name="eCGF",
    version="0.1",
    packages=find_packages() + ['cgf_pred'],
    package_dir={'cgf_pred': 'cgf_pred'},
    package_data={
        'cgf_pred': ['fastas/*.fasta', 'csvs/*.txt']
    },
    python_requires='>=3.6.0',
    entry_points={
        'console_scripts': [
            'eCGF = cgf_pred.__main__:main'
        ]
    },
)
