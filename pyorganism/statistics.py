#!/usr/bin/env python
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


__all__ = ["compute_zscore"]


import warnings
#import logging
import numpy

#from . import miscellaneous as misc
#from .errors import PyOrganismError


#LOGGER = logging.getLogger(__name__)
#LOGGER.addHandler(misc.NullHandler())


def compute_zscore(obs, random_stats):
    """
    Parameters
    ----------
    obs: numeral
        original observation
    random_stats : iterable
        same observable in randomised versions
    """
    if len(random_stats) == 0:
        return numpy.nan
    random_stats = numpy.ma.masked_invalid(random_stats)
    if random_stats.mask.all():
        warnings.warn("invalid null model values")
        return numpy.nan
    mean = numpy.mean(random_stats[numpy.logical_not(random_stats.mask)])
    std = numpy.std(random_stats[numpy.logical_not(random_stats.mask)])
    nominator = obs - mean
    if nominator == 0.0:
        return nominator
    if std == 0.0:
        if nominator < 0.0:
            return -numpy.inf
        else:
            return numpy.inf
    else:
        return (nominator / std)

