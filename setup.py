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


from os.path import join

from setuptools import (setup, Extension)
from Cython.Distutils import build_ext


if __name__ == "__main__":
    # continuous
    sources = ["_continuous.pyx", "continuous.c"]
    c_path = join("pyorganism", "regulation", "src")
    continuous = Extension("pyorganism.regulation._continuous",
        sources=[join(c_path, src) for src in sources],
        include_dirs=[c_path]
    )

    setup(
        name="pyorganism",
        version="0.1",
        description="analyse organisational principles in living organisms",
        author="Moritz Emanuel Beber",
        author_email="moritz (dot) beber (at) gmail (dot) com",
        url="http://github.com/Midnighter/pyorganism",
        packages=["pyorganism",
            "pyorganism.io",
            "pyorganism.metabolism",
            "pyorganism.regulation",
        ],
    #    package_data = {"pyorganism": ["data/*.xml", "data/*.txt", "data/*.tsv"]},
        ext_modules=[continuous],
        cmdclass={"build_ext": build_ext}
    )

