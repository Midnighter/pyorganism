# -*- coding: utf-8 -*-


"""
===========================
Regulatory Control Measures
===========================

:Authors:
    Moritz Emanuel Beber
:Date:
    2013-05-10
:Copyright:
    Copyright(c) 2013 Jacobs University of Bremen. All rights reserved.
:File:
    control.py
"""


__all__ = ["digital_control", "digital_ctc", "continuous_digital_ctc",
        "analog_control", "analog_ctc", "continuous_analog_ctc",
        "metabolic_coherence_ratio", "metabolic_coherence",
        "continuous_difference_coherence",
        "continuous_abs_coherence",
        "continuous_functional_coherence"]


import logging
import random
import re

import numpy

from copy import copy
from itertools import izip
from collections import defaultdict

from .networks import to_operon_based
from . import elements as elem
from .. import miscellaneous as misc
#from ..errors import PyOrganismError
from ..statistics import compute_zscore


LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(misc.NullHandler())

OPTIONS = misc.OptionsManager.get_instance()


def setup_trn(trn, active):
    if trn is None or trn.size() == 0:
        LOGGER.warn("empty transcriptional regulatory network")
        return numpy.nan
    active_nodes = set([gene.regulatory_product if\
        isinstance(gene.regulatory_product, elem.TranscriptionFactor) else\
        gene for gene in active])
    original = effective_network(trn, active_nodes)
    if original.size() == 0:
        LOGGER.warn("empty effective network")
        return numpy.nan
    return original

def setup_continuous_operon_based(network, elem2level):
    op2level = defaultdict(list)
    for (elem, level) in elem2level.iteritems():
        for op in elem.get_operons():
            op2level[op].append(level)
    ops = op2level.keys()
    levels = [numpy.mean(op2level[op]) for op in ops]
    op2level = dict(izip(ops, levels))
    op_net = to_operon_based(network)
    active_net = op_net.subgraph(ops)
    if len(active_net) == 0 or active_net.size() == 0:
        LOGGER.warn("empty operon network")
        return numpy.nan
    return (active_net, op2level)

def digital_control(trn, active, **kw_args):
    """
    Compute the digital control of an effective TRN using a list of
    differentially expressed genes.

    Parameters
    ----------
    trn: nx.(Multi)DiGraph
        Static transcriptional regulatory network.
    active: iterable
        An iterable with gene instances that are differentially expressed.

    Warning
    -------
    Gene instances not in the TRN are silently ignored.

    References
    ----------
    [1] Marr, C., Geertz, M., Hütt, M.-T., Muskhelishvili, G., 2008.
        Dissecting the logical types of network control in gene expression profiles.
        BMC Systems Biology 2, 18.
    """
    original = setup_trn(trn, active)
    if original is numpy.nan:
        return original
    return discrete_total_ratio(original)

def digital_ctc(trn, active, random_num=1E04, return_sample=False, **kw_args):
    """
    Compute the digital control type confidence of an effective TRN using a list
    of differentially expressed genes.

    Parameters
    ----------
    trn: nx.(Multi)DiGraph
        Static transcriptional regulatory network.
    active: iterable
        An iterable with gene instances that are differentially expressed.

    Warning
    -------
    Gene instances not in the TRN are silently ignored.

    Notes
    -----
    In the null model the number of transcription factors (nodes with
    out-degree > 0) are kept constant.

    References
    ----------
    [1] Marr, C., Geertz, M., Hütt, M.-T., Muskhelishvili, G., 2008.
        Dissecting the logical types of network control in gene expression profiles.
        BMC Systems Biology 2, 18.
    """
    random_num = int(random_num)
    original = setup_trn(trn, active)
    if original is numpy.nan:
        return original
    # new null model separates TFs and genes
    t_factors = set(node for node in trn if isinstance(node, elem.TranscriptionFactor))
    genes = set(node for node in trn if isinstance(node, elem.Gene))
    # separate numbers
    tf_num = sum(int(isinstance(item, elem.TranscriptionFactor)) for item in original)
    gene_num = len(original) - tf_num
    LOGGER.info("picked %d transcription factors", tf_num)
    samples = [trn_sample(trn, t_factors, tf_num, genes, gene_num) for i in range(random_num)]
    orig_ratio = discrete_total_ratio(original)
    z_score = compute_zscore(orig_ratio, samples)
    if return_sample:
        return (z_score, samples)
    else:
        return z_score

def continuous_digital_ctc(trn, active, expr_levels, random_num=1E04,
        return_sample=False, **kw_args):
    """
    Compute a continuous digital control type confidence of given gene
    expression levels in the effective TRN.

    Parameters
    ----------
    trn: nx.(Multi)DiGraph
        Static transcriptional regulatory network.
    active: iterable
        An iterable with gene instances.
    expr_levels: iterable
        An iterable in the same order as ``active`` with expression levels.

    Warning
    -------
    Gene instances not in the TRN are silently ignored.

    Notes
    -----
    The null model randomises expression levels but leaves the topology
    unchanged.

    References
    ----------
    [1] Marr, C., Geertz, M., Hütt, M.-T., Muskhelishvili, G., 2008.
        Dissecting the logical types of network control in gene expression profiles.
        BMC Systems Biology 2, 18.
    """
    random_num = int(random_num)
    original = setup_trn(trn, active)
    if original is numpy.nan:
        return original
    gene2level = dict(izip(active, expr_levels))
    active = original.nodes()
    orig_levels = numpy.zeros(len(active), dtype=float)
    for (i, node) in enumerate(active):
        if isinstance(node, elem.TranscriptionFactor):
            orig_levels[i] = numpy.mean([gene2level[gene] for gene in\
                    node.coded_from if gene in gene2level])
        else:
            orig_levels[i] = gene2level[node]
    orig2level = dict(izip(active, orig_levels))
    (op_net, op2level) = setup_continuous_operon_based(original, orig2level)
    # in TRN structure the out-hubs and spokes differentiation matters
    out_ops = set(node for (node, deg) in op_net.out_degree_iter() if deg > 0)
    in_ops = set(op_net.nodes_iter()).difference(out_ops)
    out_levels = [op2level[op] for op in out_ops]
    in_levels = [op2level[op] for op in in_ops]
    sample = [continuous_trn_operon_sampling(op_net, out_ops, out_levels,
            in_ops, in_levels) for i in xrange(random_num)]
    orig_ratio = continuous_functional_coherence(op_net, op2level)
    z_score = compute_zscore(orig_ratio, sample)
    if return_sample:
        return (z_score, sample)
    else:
        return z_score

def setup_gpn(gpn, active):
    if gpn is None or gpn.size() == 0:
        LOGGER.warn("empty gene proximity network")
        return numpy.nan
    original = effective_network(gpn, active)
    if original.size() == 0:
        LOGGER.warn("empty effective network")
        return numpy.nan

def analog_control(gpn, active, **kw_args):
    """
    Compute the analog control from an effective GPN.

    Parameters
    ----------
    active: iterable
        An iterable with names of genes that are active in a specific
        condition.

    Warning
    -------
    Unknown gene names are silently ignored.

    References
    ----------
    [1] Marr, C., Geertz, M., Hütt, M.-T., Muskhelishvili, G., 2008.
        Dissecting the logical types of network control in gene expression profiles.
        BMC Systems Biology 2, 18.
    """
    original = setup_gpn(gpn, active)
    if original is numpy.nan:
        return original
    return discrete_total_ratio(gpn)

def analog_ctc(gpn, active, random_num=1E04, return_sample=False, **kw_args):
    """
    Compute the analog control from an effective GPN.

    Parameters
    ----------
    active: iterable
        An iterable with names of genes that are active in a specific
        condition.

    Warning
    -------
    Unknown gene names are silently ignored.

    References
    ----------
    [1] Marr, C., Geertz, M., Hütt, M.-T., Muskhelishvili, G., 2008.
        Dissecting the logical types of network control in gene expression profiles.
        BMC Systems Biology 2, 18.
    """
    random_num = int(random_num)
    original = setup_gpn(gpn, active)
    if original is numpy.nan:
        return original
    size = len(original)
    sample = [active_sample(gpn, size) for i in xrange(random_num)]
    orig_ratio = discrete_total_ratio(original)
    z_score = compute_zscore(orig_ratio, sample)
    if return_sample:
        return (z_score, sample)
    else:
        return z_score

def continuous_analog_ctc(gpn, active, expr_levels, random_num=1E04,
        return_sample=False, **kw_args):
    """
    Compute the analog control from an effective GPN.

    Parameters
    ----------
    active: iterable
        An iterable with names of genes.
    expr_levels: iterable
        An iterable in the same order as `active` with expression levels
        (floats).

    Warning
    -------
    Unknown gene names are silently ignored.

    References
    ----------
    [1] Marr, C., Geertz, M., Hütt, M.-T., Muskhelishvili, G., 2008.
        Dissecting the logical types of network control in gene expression profiles.
        BMC Systems Biology 2, 18.
    """
    random_num = int(random_num)
    original = setup_gpn(gpn, active)
    if original is numpy.nan:
        return original
    orig2level = dict(izip(active, expr_levels))
    (op_net, op2level) = setup_continuous_operon_based(original, orig2level)
    active = op2level.keys()
    levels = [op2level[op] for op in active]
    sample = [continuous_gpn_operon_sampling(op_net, active, levels) for i in\
            xrange(random_num)]
    orig_ratio = continuous_abs_coherence(op_net, op2level)
    z_score = compute_zscore(orig_ratio, sample)
    if return_sample:
        return (z_score, sample)
    else:
        return z_score

def metabolic_coherence_ratio(metabolic_network, active, bnumber2gene, rxn_centric=None,
        **kw_args):
    """
    Compute the metabolic coherence ratio (MCR) from an effective metabolic
    network.

    Parameters
    ----------
    active: iterable
        An iterable with actively expressed genes.

    Warning
    -------
    Unknown gene names are silently ignored.

    References
    ----------
    [1] Sonnenschein, N., Geertz, M., Muskhelishvili, G., Hütt, M.-T., 2011.
        Analog regulation of metabolic demand.
        BMC Syst Biol 5, 40.

    """
    if rxn_centric is None:
        if metabolic_network is None:
            LOGGER.warn("no metabolic network")
            return numpy.nan
        else:
            rxn_centric = metabolic_network.to_reaction_centric()
    if rxn_centric.size() == 0:
        LOGGER.warn("empty metabolic network")
        return numpy.nan
    bpattern = re.compile(r"b\d{4}")
    active_reactions = list()
    # evaluate whether a reaction can be active due to gene expression
    for rxn in rxn_centric:
        info = rxn.notes["gene_association"]
        if info:
            matches = bpattern.findall(info)
            if not matches:
                continue
            for mobj in matches:
                info = info.replace(mobj, str(bnumber2gene[mobj] in active))
            activity = eval(info)
            if activity:
                active_reactions.append(rxn)
    subnet = effective_network(rxn_centric, active_reactions)
    if len(subnet) == 0:
        LOGGER.warn("empty effective network")
        return numpy.nan
    return discrete_total_ratio(subnet)

def metabolic_coherence(self, active, bnumber2gene, rxn_centric=None,
        random_num=1E04, return_sample=False, **kw_args):
    """
    Compute the metabolic coherence (MC) from an effective metabolic
    network.

    Parameters
    ----------
    active: iterable
        An iterable with actively expressed genes.

    Warning
    -------
    Unknown gene names are silently ignored.

    References
    ----------
    [1] Sonnenschein, N., Geertz, M., Muskhelishvili, G., Hütt, M.-T., 2011.
        Analog regulation of metabolic demand.
        BMC Syst Biol 5, 40.

    """
    if rxn_centric is None:
        if self.metabolic_network is None:
            LOGGER.warn("no metabolic network")
            return numpy.nan
        else:
            rxn_centric = self.metabolic_network.to_reaction_centric()
    if rxn_centric.size() == 0:
        LOGGER.warn("empty metabolic network")
        return numpy.nan
    bpattern = re.compile(r"b\d{4}")
    active_reactions = list()
    # evaluate whether a reaction can be active due to gene expression
    for rxn in rxn_centric:
        info = rxn.notes["gene_association"]
        if info:
            matches = bpattern.findall(info)
            if not matches:
                continue
            for mobj in matches:
                info = info.replace(mobj, str(bnumber2gene[mobj] in active))
            activity = eval(info)
            if activity:
                active_reactions.append(rxn)
    original = effective_network(rxn_centric, active_reactions)
    size = len(original)
    if size == 0:
        LOGGER.warn("empty effective network")
        return numpy.nan
    sample = [active_sample(rxn_centric, size) for i in xrange(int(random_num))]
    z_score = compute_zscore(discrete_total_ratio(original), sample)
    if return_sample:
        return (z_score, sample)
    else:
        return z_score

def effective_network(network, active):
    """
    Return the effective network imposed by a subset of nodes.
    """
    subnet = network.subgraph(active)
    size = len(subnet)
    LOGGER.info("{0:d}/{1:d} node(s) effective network - {2:d} entities ignored"\
            .format(size, len(network), len(active) - size))
    return subnet

def discrete_marr_ratio(network):
    control = sum(1 for (node, deg) in network.degree_iter() if deg > 0)
    if control == len(network):
        return numpy.nan
    else:
        return control / float(len(network) - control)

def discrete_total_ratio(network):
    control = sum(1 for (node, deg) in network.degree_iter() if deg > 0)
    return control / float(len(network))

def active_sample(network, size):
    """
    Sample from affected genes as null model.
    """
    sample = random.sample(network.nodes(), size)
    subnet = network.subgraph(sample)
    result = discrete_total_ratio(subnet)
    return result

def trn_sample(trn, tfs, tf_num, genes, gene_num):
    """
    Sample from affected genes and transcription factors as null model.
    """
    local_tfs = random.sample(tfs, tf_num)
    local_genes = random.sample(genes, gene_num)
    subnet = trn.subgraph(local_tfs + local_genes)
    result = discrete_total_ratio(subnet)
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
        in_levels):
    # select out- and in-ops and then compute similarity based on different
    # null-models
    op2level = dict(izip(out_ops, numpy.random.shuffle(out_levels)))
    op2level.update(izip(in_ops, numpy.random.shuffle(in_levels)))
    return continuous_functional_coherence(op_net, op2level)

def continuous_gpn_operon_sampling(op_net, active, levels):
    op2level = dict(izip(op_net, numpy.random.shuffle(levels)))
    return continuous_abs_coherence(op_net, op2level)

def continuous_difference_coherence(network, elem2level):
    return sum([elem2level[u] - elem2level[v] for (u, v) in\
            network.edges_iter()])

def continuous_abs_coherence(network, elem2level):
    return sum([1.0 - abs(elem2level[u] - elem2level[v]) for (u, v) in\
            network.edges_iter()])

def continuous_functional_coherence(network, elem2level):
    total = 0.0
    for (u, v, k) in network.edges_iter(keys=True):
        if k == 1:
            total += 1.0 - abs(elem2level[u] - elem2level[v])
        elif k == -1:
            total += abs(elem2level[u] - elem2level[v])
    return total

