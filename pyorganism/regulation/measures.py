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

from copy import copy
from itertools import izip

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
    control = sum(1 for (node, deg) in network.degree_iter() if deg > 0)
    if control == len(network):
        return numpy.nan
    else:
        return control / (len(network) - control)

def discrete_total_ratio(network):
    control = sum(1 for (node, deg) in network.degree_iter() if deg > 0)
    return control / len(network)

def active_sample(network, size, evaluate=discrete_total_ratio):
    """
    Sample from affected genes as null model.
    """
    sample = random.sample(network.nodes(), size)
    subnet = network.subgraph(sample)
    result = evaluate(subnet)
    return result

def trn_sample(trn, tfs, tf_num, genes, gene_num, evaluate=discrete_total_ratio):
    """
    Sample from affected genes and transcription factors as null model.
    """
    local_tfs = random.sample(tfs, tf_num)
    local_genes = random.sample(genes, gene_num)
    subnet = trn.subgraph(local_tfs + local_genes)
    result = evaluate(subnet)
    return result

def jack_replace(active, replacement, replace_num):
    positions = random.sample(xrange(len(active)), replace_num)
    tmp = copy(active)
    replace = random.sample(replacement, replace_num)
    for (i, k) in enumerate(positions):
        tmp[k] = replace[i]
    return tmp

def jackknife(active, remove_num):
    positions = random.sample(xrange(len(active)), remove_num)
    tmp = copy(active)
    for i in positions:
        del tmp[i]
    return tmp

def robustness(control_type, active, fraction=0.1,
        random_num=1E04, control_num=1E04, **kw_args):
    """
    Cut a fraction of active genes.
    """
    kw_args["random_num"] = control_num
    size = len(active)
    cut_num = int(round(size * fraction))
    LOGGER.info("cutting {0:d}/{1:d} genes".format(cut_num, size))
    samples = [jackknife(active, cut_num) for i in xrange(int(random_num))]
    distribution = [control_type(sample, **kw_args) for sample in samples]
    return distribution

def robustness_with_replacement(control_type, active, replacement, fraction=0.1,
        random_num=1E04, control_num=1E04, **kw_args):
    """
    Replace a fraction of the active genes with other genes.
    """
    kw_args["random_num"] = control_num
    size = len(active)
    replace_num = int(round(size * fraction))
    LOGGER.info("replacing {0:d}/{1:d} genes".format(replace_num, size))
    samples = [jack_replace(active, replacement, replace_num)\
            for i in range(int(random_num))]
    distribution = [control_type(sample, **kw_args) for sample in samples]
    return distribution

def continuous_trn_operon_sampling(op_net, out_ops, out_levels, in_ops,
        in_levels, evaluate):
    # select out- and in-ops and then compute similarity based on different
    # null-models
    numpy.random.shuffle(out_levels)
    numpy.random.shuffle(in_levels)
    op2level = dict(izip(out_ops, out_levels))
    op2level.update(izip(in_ops, in_levels))
    return evaluate(op_net, op2level)

def delayed_continuous_trn_operon_sampling(op_net, active, levels,
        delayed_levels, evaluate):
    numpy.random.shuffle(levels)
    numpy.random.shuffle(delayed_levels)
    op2level = dict(izip(active, levels))
    delayed_op2level = dict(izip(active, delayed_levels))
    return evaluate(op_net, op2level, delayed_op2level)

def continuous_gpn_operon_sampling(op_net, active, levels, evaluate):
    numpy.random.shuffle(levels)
    op2level = dict(izip(op_net, levels))
    return evaluate(op_net, op2level)

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

def delayed_continuous_difference_coherence(network, elem2level, delayed2level):
    return sum([elem2level[u] - delayed2level[v] for (u, v) in\
            network.edges_iter()]) / network.size()

def delayed_continuous_abs_coherence(network, elem2level):
    return sum([1.0 - abs(elem2level[u] - elem2level[v]) for (u, v) in\
            network.edges_iter()]) / network.size()

def delayed_continuous_functional_coherence(network, elem2level, delayed2level):
    control = 0.0
    total = 0.0
    for (u, v, k) in network.edges_iter(keys=True):
        if k == 1:
            control += 1.0 - abs(elem2level[u] - delayed2level[v])
            total += 1.0
        elif k == -1:
            control += abs(elem2level[u] - delayed2level[v])
            total += 1.0
    return control / total

def delayed_continuous_functional_comparison(network, elem2level, delayed2level):
    control = 0.0
    total = 0.0
    for (u, v, k) in network.edges_iter(keys=True):
        if k == 1:
            control += _activating(elem2level[u], delayed2level[v])
            total += 1.0
        elif k == -1:
            control += _inhibiting(elem2level[u], delayed2level[v])
            total += 1.0
    return control / total

