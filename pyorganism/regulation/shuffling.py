# -*- coding: utf-8 -*-


from __future__ import (absolute_import, unicode_literals)


"""
==============================
Regulatory Control Null Models
==============================

:Authors:
    Moritz Emanuel Beber
:Date:
    2014-10-30
:Copyright:
    Copyright(c) 2014 Jacobs University of Bremen. All rights reserved.
:File:
    shuffling.py
"""


import logging
import random

import numpy as np
from numpy.random import shuffle
from builtins import zip

from .. import miscellaneous as misc


__all__ = [
    "active_sample",
    "fixed_regulator_sample",
    "continuous_sample",
    "continuous_fixed_regulator_sample",
    "delayed_continuous_sample",
    "timeline_sample"
]


LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(misc.NullHandler())

OPTIONS = misc.OptionsManager.get_instance()


def active_sample(network, size, evaluate):
    """
    Sample from affected genes as null model.
    """
    sample = random.sample(network.nodes(), size)
    subnet = network.subgraph(sample)
    result = evaluate(subnet)
    return result

def fixed_regulator_sample(reference, regulators, reg_num, slaves, slave_num,
        evaluate):
    """
    Sample from affected genes and transcription factors as null model.
    """
    eff_regs = random.sample(regulators, reg_num)
    eff_slaves = random.sample(slaves, slave_num)
    subnet = reference.subgraph(eff_regs + eff_slaves)
    return evaluate(subnet)

def continuous_sample(net, active, levels, evaluate):
    shuffle(levels)
    node2level = {node: lvl for (node, lvl) in zip(active, levels)}
    return evaluate(net, node2level)

def continuous_fixed_regulator_sample(net, regulators, reg_levels, slaves,
        slave_levels, evaluate):
    # select out- and in-ops and then compute similarity based on different
    # null-models
    shuffle(reg_levels)
    shuffle(slave_levels)
    node2level = {node: lvl for (node, lvl) in zip(regulators, reg_levels)}
    node2level.update(zip(slaves, slave_levels))
    return evaluate(net, node2level)

def delayed_continuous_sample(net, active, levels, delayed_levels, evaluate):
    shuffle(levels)
    shuffle(delayed_levels)
    node2level = {node: lvl for (node, lvl) in zip(active, levels)}
    node2delayed = {node: lvl for (node, lvl) in zip(active, delayed_levels)}
    return evaluate(net, node2level, node2delayed)

def timeline_sample(series, num):
    rnd_cols = np.arange(series.shape[1])
    for i in range(num):
        shuffle(rnd_cols)
        yield np.ascontiguousarray(series[:, rnd_cols])

