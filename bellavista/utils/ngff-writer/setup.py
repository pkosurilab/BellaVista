#!/usr/bin/env python

from distutils.core import setup

packages = ["ngff_writer"]

package_data = {"": ["*"]}

install_requires = [
    "dask>=2.11.0",
    "dask_image",
    "numpy",
    "scikit-image",
    "zarr",
]

extras_require = {"test": ["pytest~=6.2"]}

setup(
    name="ngff-writer",
    version="0.0.1",
    description="Library for writing and reading OME-NGFF images",
    author="Alexandrov Team, EMBL",
    author_email=None,
    url=None,
    packages=packages,
    package_data=package_data,
    install_requires=install_requires,
    extras_require=extras_require,
    python_requires=">=3.8",
)
