from setuptools import setup, find_packages

setup(
        name="jwst-destripe-driver",
        version="0.0",
        packages=find_packages(),
        install_requires=[
                'numpy',
                'argparse'
        ]
)  
