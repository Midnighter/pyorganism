#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
====================
Micro Array Data I/O
====================

:Authors:
    Moritz Emanuel Beber
:Date:
    2014-01-14
:Copyright:
    Copyright(c) 2014 Jacobs University of Bremen. All rights reserved.
:File:
    microarray.py
"""


__all__ = ["read_microarray", "read_interpolated"]


import logging

import numpy as np
import pandas

from .. import miscellaneous as misc


LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(misc.NullHandler())


def read_microarray(filename):
    return pandas.read_table(filename, sep="\t", names=["name", "blattner",
            "ratio (A/B)", "p-value", "function"], skiprows=1)

def read_interpolated(filename, mutants=[]):
    df = pandas.read_table(filename, sep="\t", index_col=0, header=0)
    df.sort_index(inplace=True)
    for ind in mutants:
        df.loc[ind] = np.nan
    return df

