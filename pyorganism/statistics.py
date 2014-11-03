# -*- coding: utf-8 -*-


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
#from .errors import PyOrganismError


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
    Normalise a np.array to values between zero and unity.

    Warning
    -------
    Beware of infinity values.
    """
    mask = np.isfinite(vec)
    if len(vec[mask]) == 0:
        return vec
    min_num = vec[mask].min()
    max_num = vec[mask].max()
    return (vec - min_num) / (max_num - min_num)

