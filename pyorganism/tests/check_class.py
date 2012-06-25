#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
========================
Metabolic Elements Tests
========================

:Authors:
    Moritz Emanuel Beber
:Date:
    2011-06-25
:Copyright:
    Copyright(c) 2012 Jacobs University of Bremen. All rights reserved.
:File:
    check_class.py
"""


import itertools

import nose.tools as nt


def check__str__(instances):
    for obj in instances:
        nt.assert_equal(str(obj.name), str(obj))

def check__unicode__(instances):
    for obj in instances:
        nt.assert_equal(unicode(obj.name), unicode(obj))

def check__repr__(instances):
    for obj in instances:
        nt.assert_equal(u"<{0}.{1} {2:d}>".format(obj.__module__,
                obj.__class__.__name__, id(obj)), repr(obj))

def check_reversible(instances, values):
    nt.assert_equal(len(instances), len(values))
    for (obj, value) in itertools.izip(instances, values):
        nt.assert_equal(obj.reversible, value)

