#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
=====================
PyOrganism SBML Tests
=====================

:Authors:
    Moritz Emanuel Beber
:Date:
    2013-09-05
:Copyright:
    Copyright(c) 2013 Jacobs University of Bremen. All rights reserved.
:File:
    test_sbml.py
"""


#import nose.tools as nt

from glob import glob

from pyorganism.io.sbml import SBMLParser


def test_parser():
    parser = SBMLParser()
    for filename in glob("models/*.xml"):
        met = parser(filename)
        print len(met.compartments)
        print len(met.compounds)
        print len(met.reactions)

