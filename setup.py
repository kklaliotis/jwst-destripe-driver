from setuptools import setup, find_packages

setup(
        name="jwst-destripe-driver",
        version="0.1.0",
        package_dir={"": "src"},
        packages=find_packages(where="src"),
        install_requires=[
                "numpy",
                "astropy",
                "pyimcom",
        ],
        entry_points={
                "console_scripts": [
                        "jwst-destripe=jwst_destripe_driver.run:main",
                ]
        },
)