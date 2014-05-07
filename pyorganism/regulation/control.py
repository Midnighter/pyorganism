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


__all__ = ["digital_control", "digital_ctc", "continuous_digital_control",
        "continuous_digital_ctc_fixed_tf_num", "delayed_continuous_digital_ctc",
        "analog_control", "analog_ctc", "continuous_analog_control",
        "continuous_analog_ctc",
        "metabolic_coherence_ratio", "metabolic_coherence"]


import logging
import re

import numpy

from itertools import izip

from . import measures as ms
from . import networks as nets
from .. import miscellaneous as misc
#from ..errors import PyOrganismError
from ..statistics import compute_zscore


LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(misc.NullHandler())

OPTIONS = misc.OptionsManager.get_instance()


def digital_control(trn, active, measure=ms.discrete_total_ratio):
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
    (original, _, _) = nets.setup_trn(trn, active)
    if numpy.isnan(original):
        return original
    return measure(original)

def digital_ctc(trn, active, random_num=1E04, return_sample=False,
        measure=ms.discrete_total_ratio):
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
    (original, _, _) = nets.setup_trn(trn, active)
    if numpy.isnan(original):
        return original
    size = len(original)
    sample = [ms.active_sample(trn, size, evaluate=measure) for i in xrange(random_num)]
    z_score = compute_zscore(measure(original), sample)
    if return_sample:
        return (z_score, sample)
    else:
        return z_score

def digital_ctc_fixed_tf_num(trn, active, random_num=1E04, return_sample=False,
        measure=ms.discrete_total_ratio):
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
    (original, t_factors, regulated) = nets.setup_trn(trn, active)
    if numpy.isnan(original):
        return original
    # separate numbers
    tf_num = len(t_factors)
    gene_num = len(regulated)
    LOGGER.info("picked %d transcription factors", tf_num)
    sample = [ms.trn_sample(trn, t_factors, tf_num, regulated, gene_num,
            evaluate=measure) for i in range(random_num)]
    z_score = compute_zscore(measure(original), sample)
    if return_sample:
        return (z_score, sample)
    else:
        return z_score

def continuous_digital_control(trn, active, expr_levels,
        measure=ms.continuous_functional_coherence):
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
    (original, t_factors, regulated) = nets.setup_trn(trn, active)
    if numpy.isnan(original):
        return original
    gene2level = dict(izip(active, expr_levels))
    orig2level = nets.setup_trn_levels(t_factors, regulated, gene2level)
    (op_net, op2level) = nets.setup_continuous_operon_based(original, orig2level)
    return measure(op_net, op2level)

def continuous_digital_ctc_fixed_tf_num(trn, active, expr_levels, random_num=1E04,
        return_sample=False, measure=ms.continuous_functional_coherence):
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
    (original, t_factors, regulated) = nets.setup_trn(trn, active)
    if numpy.isnan(original):
        return original
    gene2level = dict(izip(active, expr_levels))
    orig2level = nets.setup_trn_levels(t_factors, regulated, gene2level)
    (op_net, op2level) = nets.setup_continuous_operon_based(original, orig2level)
    # in TRN structure the out-hubs and spokes differentiation matters
    out_ops = set(node for (node, deg) in op_net.out_degree_iter() if deg > 0)
    if len(out_ops) == 0:
        LOGGER.error("no out-hubs in operon-based TRN")
    in_ops = set(op_net.nodes_iter()).difference(out_ops)
    if len(in_ops) == 0:
        LOGGER.error("no targets in operon-based TRN")
    out_levels = [op2level[op] for op in out_ops]
    if len(out_levels) == 0:
        LOGGER.error("no out-hub expression levels")
    in_levels = [op2level[op] for op in in_ops]
    if len(in_ops) == 0:
        LOGGER.error("no target expression levels")
    sample = [ms.continuous_trn_operon_sampling(op_net, out_ops, out_levels,
            in_ops, in_levels, evaluater=measure) for i in xrange(random_num)]
    orig_ratio = measure(op_net, op2level)
    z_score = compute_zscore(orig_ratio, sample)
    if return_sample:
        return (z_score, sample)
    else:
        return z_score

def delayed_continuous_digital_ctc(trn, active, expr_levels,
        delayed_expr_levels, random_num=1E04, return_sample=False,
        measure=ms.continuous_functional_coherence):
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
    (original, t_factors, regulated) = nets.setup_trn(trn, active)
    if numpy.isnan(original):
        return original
    gene2level = dict(izip(active, expr_levels))
    orig2level = nets.setup_trn_levels(t_factors, regulated, gene2level)
    gene2level = dict(izip(active, delayed_expr_levels))
    delayed2level = nets.setup_trn_levels(t_factors, regulated, gene2level)
    (op_net, op2level) = nets.setup_continuous_operon_based(original, orig2level)
    (op_net, delayed_op2level) = nets.setup_continuous_operon_based(original,
            delayed2level)
    active_nodes = op_net.nodes()
# Ignoring out-hubs for time delayed analysis for now
#    # in TRN structure the out-hubs and spokes differentiation matters
#    out_ops = set(node for (node, deg) in op_net.out_degree_iter() if deg > 0)
#    if len(out_ops) == 0:
#        LOGGER.error("no out-hubs in operon-based TRN")
#    in_ops = set(op_net.nodes_iter()).difference(out_ops)
#    if len(in_ops) == 0:
#        LOGGER.error("no targets in operon-based TRN")
#    out_levels = [op2level[op] for op in out_ops]
#    if len(out_levels) == 0:
#        LOGGER.error("no out-hub expression levels")
#    in_levels = [op2level[op] for op in in_ops]
#    if len(in_ops) == 0:
#        LOGGER.error("no target expression levels")
    levels = [op2level[node] for node in active_nodes]
    delayed_levels = [delayed_op2level[node] for node in active_nodes]
    sample = [ms.delayed_continuous_trn_operon_sampling(op_net, active_nodes,
            levels, delayed_levels, measure) for i in xrange(random_num)]
    orig_ratio = measure(op_net, op2level, delayed_op2level)
    z_score = compute_zscore(orig_ratio, sample)
    if return_sample:
        return (z_score, sample)
    else:
        return z_score

def analog_control(gpn, active):
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
    original = nets.setup_gpn(gpn, active)
    if numpy.isnan(original):
        return original
    return ms.discrete_total_ratio(gpn)

def analog_ctc(gpn, active, random_num=1E04, return_sample=False,
        measure=ms.discrete_total_ratio):
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
    original = nets.setup_gpn(gpn, active)
    if numpy.isnan(original):
        return original
    size = len(original)
    sample = [ms.active_sample(gpn, size, evaluate=measure) for i in xrange(random_num)]
    z_score = compute_zscore(measure(original), sample)
    if return_sample:
        return (z_score, sample)
    else:
        return z_score

def continuous_analog_control(gpn, active, expr_levels,
        measure=ms.continuous_abs_coherence):
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
    original = nets.setup_gpn(gpn, active)
    if numpy.isnan(original):
        return original
    gene2level = dict(izip(active, expr_levels))
    (op_net, op2level) = nets.setup_continuous_operon_based(original, gene2level)
    return measure(op_net, op2level)

def continuous_analog_ctc(gpn, active, expr_levels, random_num=1E04,
        return_sample=False, measure=ms.continuous_abs_coherence):
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
    original = nets.setup_gpn(gpn, active)
    if numpy.isnan(original):
        return original
    gene2level = dict(izip(active, expr_levels))
    (op_net, op2level) = nets.setup_continuous_operon_based(original, gene2level)
    active = op_net.nodes()
    levels = [op2level[op] for op in active]
    sample = [ms.continuous_gpn_operon_sampling(op_net, active, levels,
        measure) for i in xrange(random_num)]
    orig_ratio = measure(op_net, op2level)
    z_score = compute_zscore(orig_ratio, sample)
    if return_sample:
        return (z_score, sample)
    else:
        return z_score

def metabolic_coherence_ratio(metabolic_network, active, bnumber2gene,
        rxn_centric=None, measure=ms.discrete_total_ratio):
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
    rxn_centric = nets.setup_metabolic(metabolic_network, rxn_centric)
    if numpy.isnan(rxn_centric):
        return rxn_centric
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
    original = ms.effective_network(rxn_centric, active_reactions)
    if len(original) == 0:
        LOGGER.warn("empty effective network")
        return numpy.nan
    return measure(original)

def metabolic_coherence(metabolic_network, active, bnumber2gene, rxn_centric=None,
        random_num=1E04, return_sample=False, measure=ms.discrete_total_ratio):
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
    random_num = int(random_num)
    rxn_centric = nets.setup_metabolic(metabolic_network, rxn_centric)
    if numpy.isnan(rxn_centric):
        return rxn_centric
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
    original = ms.effective_network(rxn_centric, active_reactions)
    size = len(original)
    if size == 0:
        LOGGER.warn("empty effective network")
        return numpy.nan
    sample = [ms.active_sample(rxn_centric, size, measure) for i in xrange(int(random_num))]
    z_score = compute_zscore(measure(original), sample)
    if return_sample:
        return (z_score, sample)
    else:
        return z_score

