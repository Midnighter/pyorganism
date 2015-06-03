#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
==================
PyOrganism Package
==================

:Authors:
    Moritz Emanuel Beber
:Date:
    2012-05-22
:Copyright:
    Copyright(c) 2012 Jacobs University of Bremen. All rights reserved.
:File:
    setup.py
"""


import sys
from os.path import join

from setuptools import (setup, Extension)
try:
    from Cython.Distutils import build_ext
except ImportError as err:
    sys.exit("Apologies, you need 'Cython' to install 'pyorganism'.")


if __name__ == "__main__":
    # continuous
    sources = ["continuous_wrapper.pyx", "continuous.c"]
    c_path = join("pyorganism", "regulation", "src")
    continuous = Extension("pyorganism.regulation.continuous_wrapper",
        sources=[join(c_path, src) for src in sources],
        include_dirs=[c_path]
    )

    setup(
        name="pyorganism",
        version="0.2.4",
        license="BSD",
        description="analyze organisational principles in living organisms",
        author="Moritz Emanuel Beber",
        author_email="moritz (dot) beber (at) gmail (dot) com",
        url="http://github.com/Midnighter/pyorganism",
        zip_safe=False,
        install_requires=[
            "future",
            "networkx",
            "numpy",
            "pandas"
        ],
        packages=["pyorganism",
            "pyorganism.io",
            "pyorganism.metabolism",
            "pyorganism.regulation"m ,
        ],
    #    package_data = {"pyorganism": ["data/*.xml", "data/*.txt", "data/*.tsv"]},
        ext_modules=[continuous],
        cmdclass={"build_ext": build_ext},
        classifiers=[
            # complete classifier list: http://pypi.python.org/pypi?%3Aaction=list_classifiers
            "Development Status :: 4 - Beta",
            "Intended Audience :: Developers",
            "License :: OSI Approved :: BSD License",
            "Natural Language :: English",
            "Operating System :: Unix",
            "Operating System :: POSIX",
            "Operating System :: Microsoft :: Windows",
            "Programming Language :: Python",
            "Programming Language :: Python :: 2.6",
            "Programming Language :: Python :: 2.7",
            "Programming Language :: Python :: 3",
            "Programming Language :: Python :: 3.3",
            "Programming Language :: Python :: 3.4",
            "Programming Language :: Python :: Implementation :: CPython",
            "Topic :: Scientific/Engineering :: Bio-Informatics",
        ],
    )

