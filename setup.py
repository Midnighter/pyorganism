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


from distutils.core import setup


setup(
    name = "pyorganism",
    version = "0.0",
    description = "analyse organisational principles in living organisms",
    author = "Moritz Emanuel Beber",
    author_email = "moritz (dot) beber (at) googlemail (dot) com",
    url = "http://github.com/Midnighter/pymetabolism",
    packages = ["pyorganism",
            "pyorganism.tests"],
    package_data = {"pyorganism": ["data/*.xml", "data/*.txt", "data/*.tsv"]},
    )

