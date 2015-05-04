# -*- coding: utf-8 -*-


from __future__ import (absolute_import, unicode_literals, division)


"""
=====================
PyOrganism Statistics
=====================

:Authors:
    Moritz Emanuel Beber
:Date:
    2012-06-01
:Copyright:
    Copyright(c) 2012 Jacobs University of Bremen. All rights reserved.
:File:
    statistics.py
"""


__all__ = ["compute_zscore", "norm_zero2unity"]


import logging

import numpy as np

from . import miscellaneous as misc


LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(misc.NullHandler())


def compute_zscore(obs, random_stats):
    """
    Parameters
    ----------
    obs: numeral
        original observation
    random_stats : iterable
        same observable in randomised versions
    """
    random_stats = np.asarray(random_stats)
    mask = np.isfinite(random_stats)
    if len(random_stats[mask]) == 0:
        LOGGER.warn("invalid null model values")
        return np.nan
    values = random_stats[mask]
    mean = np.mean(values)
    std = np.std(values)
    nominator = obs - mean
    if nominator == 0.0:
        return nominator
    if std == 0.0:
        if nominator < 0.0:
            return -np.inf
        else:
            return np.inf
    else:
        return (nominator / std)

def norm_zero2unity(vec):
    """
    Normalise a np.array to values between zero and unity. Except when all
    values are equal or NaN.

    Warning
    -------
    Beware of infinity values.
    """
    mask = np.isfinite(vec)
    if len(vec[mask]) == 0:
        return vec
    elif (vec == 0).all():
        return vec
    min_val = vec[mask].min()
    max_val = vec[mask].max()
    if min_val == max_val:
        return np.zeros_like(vec)
    res = vec - min_val
    res /= max_val - min_val
    return res

