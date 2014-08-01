# -*- coding: utf-8 -*-


"""
===========================
Regulatory Control Measures
===========================

:Authors:
    Moritz Emanuel Beber
:Date:
    2013-10-11
:Copyright:
    Copyright(c) 2013 Jacobs University of Bremen. All rights reserved.
:File:
    measures.py
"""


from __future__ import division


import logging
import random

import numpy

from itertools import izip

from numpy.random import shuffle

from .. import miscellaneous as misc
#from ..errors import PyOrganismError


__all__ = ["discrete_marr_ratio",
        "discrete_total_ratio",
        "continuous_abs_coherence",
        "continuous_difference_coherence",
        "continuous_abs_difference_coherence",
        "continuous_functional_coherence",
        "continuous_functional_comparison",
        "delayed_continuous_difference_coherence",
        "delayed_continuous_abs_coherence",
        "delayed_continuous_functional_coherence",
        "delayed_continuous_functional_comparison"]


LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(misc.NullHandler())

OPTIONS = misc.OptionsManager.get_instance()


def discrete_marr_ratio(network):
    control = sum(deg > 0 for (node, deg) in network.degree_iter())
    if control == len(network):
        return numpy.nan
    else:
        return control / (len(network) - control)

def discrete_total_ratio(network):
    control = sum(deg > 0 for (node, deg) in network.degree_iter())
    return control / len(network)

def active_sample(network, size, evaluate=discrete_total_ratio):
    """
    Sample from affected genes as null model.
    """
    sample = random.sample(network.nodes(), size)
    subnet = network.subgraph(sample)
    result = evaluate(subnet)
    return result

def fixed_regulator_sample(reference, regulators, reg_num, slaves, slave_num,
        evaluate=discrete_total_ratio):
    """
    Sample from affected genes and transcription factors as null model.
    """
    eff_regs = random.sample(regulators, reg_num)
    eff_slaves = random.sample(slaves, slave_num)
    subnet = reference.subgraph(eff_regs + eff_slaves)
    return evaluate(subnet)

def continuous_sample(net, active, levels, evaluate):
    shuffle(levels)
    node2level = {node: lvl for (node, lvl) in izip(active, levels)}
    return evaluate(net, node2level)

def continuous_fixed_regulator_sample(net, regulators, reg_levels, slaves,
        slave_levels, evaluate):
    # select out- and in-ops and then compute similarity based on different
    # null-models
    shuffle(reg_levels)
    shuffle(slave_levels)
    node2level = {node: lvl for (node, lvl) in izip(regulators, reg_levels)}
    node2level.update(izip(slaves, slave_levels))
    return evaluate(net, node2level)

def delayed_continuous_sample(net, active, levels, delayed_levels, evaluate):
    shuffle(levels)
    shuffle(delayed_levels)
    node2level = {node: lvl for (node, lvl) in izip(active, levels)}
    node2delayed = {node: lvl for (node, lvl) in izip(active, delayed_levels)}
    return evaluate(net, node2level, node2delayed)

def continuous_abs_coherence(network, elem2level):
    return sum([1.0 - abs(elem2level[u] - elem2level[v]) for (u, v) in\
            network.edges_iter()]) / network.size()

def continuous_difference_coherence(network, elem2level):
    return sum([elem2level[u] - elem2level[v] for (u, v) in\
            network.edges_iter()]) / network.size()

def continuous_abs_difference_coherence(network, elem2level):
    return sum([abs(elem2level[u] - elem2level[v]) for (u, v) in\
            network.edges_iter()]) / network.size()

def continuous_functional_coherence(network, elem2level):
    control = 0.0
    total = 0.0
    for (u, v, k) in network.edges_iter(keys=True):
        if k == 1:
            control += 1.0 - abs(elem2level[u] - elem2level[v])
            total += 1.0
        elif k == -1:
            control += abs(elem2level[u] - elem2level[v])
            total += 1.0
    return control / total

def _activating(rate_a, rate_b):
    if (rate_a > 0):
        return float(rate_b >= 0)
    elif (rate_a <= 0): # could use else here
        return float(rate_b <= 0)

def _inhibiting(rate_a, rate_b):
    if (rate_a > 0):
        return float(rate_b <= 0)
    elif (rate_a <= 0): # could use else here
        return float(rate_b >= 0)

def continuous_functional_comparison(network, elem2level):
    control = 0.0
    total = 0.0
    for (u, v, k) in network.edges_iter(keys=True):
        if k == 1:
            control += _activating(elem2level[u], elem2level[v])
            total += 1.0
        elif k == -1:
            control += _inhibiting(elem2level[u], elem2level[v])
            total += 1.0
    return control / total

def delayed_continuous_difference_coherence(network, elem2level, elem2delayed):
    return sum([elem2level[u] - elem2delayed[v] for (u, v) in\
            network.edges_iter()]) / network.size()

def delayed_continuous_abs_coherence(network, elem2level, elem2delayed):
    return sum([1.0 - abs(elem2level[u] - elem2delayed[v]) for (u, v) in\
            network.edges_iter()]) / network.size()

def delayed_continuous_functional_coherence(network, elem2level, elem2delayed):
    control = 0.0
    total = 0.0
    for (u, v, k) in network.edges_iter(keys=True):
        if k == 1:
            control += 1.0 - abs(elem2level[u] - elem2delayed[v])
            total += 1.0
        elif k == -1:
            control += abs(elem2level[u] - elem2delayed[v])
            total += 1.0
    return control / total

def delayed_continuous_functional_comparison(network, elem2level, elem2delayed):
    control = 0.0
    total = 0.0
    for (u, v, k) in network.edges_iter(keys=True):
        if k == 1:
            control += _activating(elem2level[u], elem2delayed[v])
            total += 1.0
        elif k == -1:
            control += _inhibiting(elem2level[u], elem2delayed[v])
            total += 1.0
    return control / total

