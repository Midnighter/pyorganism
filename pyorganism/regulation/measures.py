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

import numpy as np

from itertools import izip
from collections import defaultdict

from .. import miscellaneous as misc
from .elements import TranscriptionFactor


__all__ = [
    "active_genes_and_tf",
    "gene_and_tf_levels",
    "active_tu",
    "tu_levels",
    "active_operons",
    "operon_levels",
    "discrete_marr_ratio",
    "discrete_total_ratio",
]


LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(misc.NullHandler())

OPTIONS = misc.OptionsManager.get_instance()


def active_genes_and_tf(active):
    """
    From differentially expressed genes, return the genes and transcription
    factors.
    """
    active = list(active)
    t_factors = list({gene.regulatory_product for gene in active if\
            isinstance(gene.regulatory_product, TranscriptionFactor)})
    return t_factors + active

def gene_and_tf_levels(active, levels):
    """
    From genes and their respective expression levels, return the genes and transcription
    factors with their mean levels.
    """
    active = list(active)
    t_factors = list({gene.regulatory_product for gene in active if\
            isinstance(gene.regulatory_product, TranscriptionFactor)})
    orig_levels = np.zeros(len(t_factors) + len(active), dtype=float)
    gene2level = dict(izip(active, levels))
    for (i, tf) in enumerate(t_factors):
        orig_levels[i] = np.mean([gene2level[gene] for gene in\
                tf.coded_from if gene in gene2level])
    for (i, gene) in enumerate(active, start=len(t_factors)):
        orig_levels[i] = gene2level[gene]
    return (t_factors + active, orig_levels)

def active_tu(active):
    """
    From differentially expressed genes, return the transcription units.
    """
    active = list(active)
    t_units = list({tu for gene in active for tu in gene.transcription_units})
    return t_units

def tu_levels(active, levels):
    """
    From genes and their respective expression levels, return the transcription
    units with their mean levels.
    """
    active = list(active)
    gene2level = dict(izip(active, levels))
    tu2level = defaultdict(list)
    for gene in active:
        # genes without associated TU are their own TU
        if len(gene.transcription_units) == 0:
            tu2level[gene].append(gene2level[gene])
        else:
            for tu in gene.transcription_units:
                tu2level[tu].append(gene2level[gene])
    t_units = tu2level.keys()
    levels = [np.mean(tu2level[tu]) for tu in t_units]
    return (t_units, levels)

def active_operons(active):
    """
    From differentially expressed genes, return the operons.
    """
    active = list(active)
    operons = list({tu for gene in active for tu in gene.operons})
    return operons

def operon_levels(active, levels):
    """
    From genes and their respective expression levels, return the operons
    with their mean levels.
    """
    active = list(active)
    gene2level = dict(izip(active, levels))
    op2level = defaultdict(list)
    for gene in active:
        # genes without associated operons are their own operon
        if len(gene.operons) == 0:
            op2level[gene].append(gene2level[gene])
        else:
            for op in gene.operons:
                op2level[op].append(gene2level[gene])
    operons = op2level.keys()
    levels = [np.mean(op2level[op]) for op in operons]
    return (operons, levels)

def discrete_marr_ratio(network):
    control = sum(deg > 0 for (node, deg) in network.degree_iter())
    if control == len(network):
        return np.nan
    else:
        return control / (len(network) - control)

def discrete_total_ratio(network):
    control = sum(deg > 0 for (node, deg) in network.degree_iter())
    return control / len(network)

