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
    test_base.py
"""


import nose.tools as nt

from . import check_class as check
from ..base import UniqueBase


def test_components():
    foo = UniqueBase("foo")
    nt.assert_equal(foo, UniqueBase("foo"))
    args = tuple([None] * 10)
    bar = UniqueBase("bar", *args)
    nt.assert_equal(bar, UniqueBase("bar"))
    kw_args = {"a": 1, "b": 2, "c": 3}
    snafu = UniqueBase("snafu", *args, **kw_args)
    nt.assert_equal(snafu, UniqueBase("snafu"))
    objs = [foo, bar, snafu, UniqueBase()]
    check.check__str__(objs)
    check.check__unicode__(objs)
    check.check__repr__(objs)

