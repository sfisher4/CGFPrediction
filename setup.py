from setuptools import setup, find_packages

setup(
    name="eCGF",
    version="0.1",
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'my_project = my_project.__main__:main'
        ]
    },
)