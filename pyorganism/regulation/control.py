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


__all__ = ["digital_control", "digital_ctc", "digital_ctc_fixed_regulators",
#        "continuous_digital_control", "continuous_digital_ctc",
#        "continuous_digital_ctc_fixed_regulators",
#        "delayed_continuous_digital_control",
#        "delayed_continuous_digital_ctc",
        "analog_control", "analog_ctc",
#        "continuous_analog_control", "continuous_analog_ctc",
        "metabolic_coherence_ratio", "metabolic_coherence"]


import logging
import re

import numpy as np

#from itertools import izip

from . import measures as ms
from . import shuffling as shuff
from . import networks as nets
from .. import miscellaneous as misc
#from ..errors import PyOrganismError
from ..statistics import compute_zscore


LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(misc.NullHandler())

OPTIONS = misc.OptionsManager.get_instance()


def digital_control(effective, measure=ms.discrete_total_ratio):
    """
    Compute the digital control of an effective transcriptional regulatory
    network (TRN).

    Parameters
    ----------
    effective: TRN or nx.(Multi)DiGraph
        Effective TRN.
    measure: callable (optional)
        Takes the effective network as its only argument and returns the
        magnitude of control within it.

    References
    ----------
    [1] Marr, C., Geertz, M., Hütt, M.-T., Muskhelishvili, G., 2008.
        Dissecting the logical types of network control in gene expression profiles.
        BMC Systems Biology 2, 18.
    """
    if effective is None or effective.size() == 0:
        return np.nan
    return measure(effective)

def digital_ctc(effective, reference, measure=ms.discrete_total_ratio,
        random_num=1E04, return_sample=False):
    """
    Compute the digital control type confidence (CTC) of an effective
    transcriptional regulatory network (TRN).

    This function computes a Z-score of digital control using a reference
    (complete) TRN as a simple null model.

    Parameters
    ----------
    effective: TRN or nx.(Multi)DiGraph
        Effective TRN.
    reference: TRN or nx.(Multi)DiGraph
        Complete TRN.
    measure: callable (optional)
        Takes the effective network as its only argument and returns the
        magnitude of control within it.
    random_num: int (optional)
        Size of the sample distribution of digital control in the null model.
    return_sample: bool (optional)
        Whether or not to return the sample distribution.

    References
    ----------
    [1] Marr, C., Geertz, M., Hütt, M.-T., Muskhelishvili, G., 2008.
        Dissecting the logical types of network control in gene expression profiles.
        BMC Systems Biology 2, 18.
    """
    random_num = int(random_num)
    if effective is None or effective.size() == 0:
        return np.nan
    size = len(effective)
    sample = [shuff.active_sample(reference, size, evaluate=measure) for i in xrange(random_num)]
    z_score = compute_zscore(measure(effective), sample)
    if return_sample:
        return (z_score, sample)
    else:
        return z_score

def digital_ctc_fixed_regulators(effective, reference, measure=ms.discrete_total_ratio,
        random_num=1E04, return_sample=False):
    """
    Compute the digital control type confidence (CTC) of an effective
    transcriptional regulatory network (TRN).

    This function computes a Z-score of digital control using an improved null
    model where the number of regulating nodes (out-degree > 0) is kept constant.

    Parameters
    ----------
    effective: TRN or nx.(Multi)DiGraph
        Effective TRN.
    reference: TRN or nx.(Multi)DiGraph
        Complete TRN.
    measure: callable (optional)
        Takes the effective network as its only argument and returns the
        magnitude of control within it.
    random_num: int (optional)
        Size of the sample distribution of digital control in the null model.
    return_sample: bool (optional)
        Whether or not to return the sample distribution.

    References
    ----------
    [1] 
    """
    random_num = int(random_num)
    if effective is None or effective.size() == 0:
        return np.nan
    (eff_regs, eff_slaves) = nets.split_regulators(effective)
    (regulators, slaves) = nets.split_regulators(reference)
    LOGGER.info("picked %d regulators", len(eff_regs))
    sample = [shuff.fixed_regulator_sample(reference, regulators, len(eff_regs),
            slaves, len(eff_slaves), evaluate=measure) for i in range(random_num)]
    z_score = compute_zscore(measure(effective), sample)
    if return_sample:
        return (z_score, sample)
    else:
        return z_score

#def continuous_digital_control(trn, active, levels,
#        measure=ms.continuous_functional_coherence):
#    """
#    Compute the continuous digital control of a transcriptional regulatory
#    network (TRN).
#
#    Uses expression levels of nodes in a TRN in order to evaluate the magnitude of
#    control.
#
#    Parameters
#    ----------
#    trn: TRN or nx.(Multi)DiGraph
#        Effective TRN.
#    active: list
#        Ordered collection of active nodes.
#    levels: list
#        Corresponding expression levels of active nodes (same shape and ordering
#        expected).
#    measure: callable (optional)
#        Takes the effective network and expression level map and returns the
#        magnitude of control.
#
#    References
#    ----------
#    [1] 
#    """
#    if trn is None or trn.size() == 0:
#        return np.nan
#    node2level = {node: lvl for (node, lvl) in izip(active, levels)}
#    return measure(trn, node2level)
#
#def continuous_digital_ctc(trn, active, levels,
#        measure=ms.continuous_functional_coherence, random_num=1E04,
#        return_sample=False):
#    """
#    Compute the continuous digital control type confidence (CTC) in a
#    transcriptional regulatory network (TRN).
#
#    This function computes a Z-score of continuous digital control using
#    expression levels of nodes in a TRN as a simple null model.
#
#    Parameters
#    ----------
#    trn: TRN or nx.(Multi)DiGraph
#        Effective TRN.
#    active: list
#        Ordered collection of active nodes.
#    levels: list
#        Corresponding expression levels of active nodes (same shape and ordering
#        expected).
#    measure: callable (optional)
#        Takes the effective network and expression level map and returns the
#        magnitude of control.
#    random_num: int (optional)
#        Size of the sample distribution of continuous digital control in the
#        null model.
#    return_sample: bool (optional)
#        Whether or not to return the sample distribution.
#
#    References
#    ----------
#    [1] 
#    """
#    random_num = int(random_num)
#    if trn is None or trn.size() == 0:
#        return np.nan
#    node2level = {node: lvl for (node, lvl) in izip(active, levels)}
#    sample = [shuff.continuous_sample(trn, active, levels, evaluate=measure)\
#            for i in xrange(random_num)]
#    z_score = compute_zscore(measure(trn, node2level), sample)
#    if return_sample:
#        return (z_score, sample)
#    else:
#        return z_score
#
#def continuous_digital_ctc_fixed_regulators(trn, active, levels, random_num=1E04,
#        return_sample=False, measure=ms.continuous_functional_coherence):
#    """
#    Compute the continuous digital control type confidence (CTC) in a
#    transcriptional regulatory network (TRN).
#
#    This function computes a Z-score of continuous digital control using
#    expression levels of nodes in a TRN and considering expression levels of
#    regulating nodes (out-degree > 0) and regulated nodes (out-degree = 0)
#    separately.
#
#    Parameters
#    ----------
#    trn: TRN or nx.(Multi)DiGraph
#        Effective TRN.
#    active: list
#        Ordered collection of active nodes.
#    levels: list
#        Corresponding expression levels of active nodes (same shape and ordering
#        expected).
#    measure: callable (optional)
#        Takes the effective network and expression level map and returns the
#        magnitude of control.
#    random_num: int (optional)
#        Size of the sample distribution of continuous digital control in the
#        null model.
#    return_sample: bool (optional)
#        Whether or not to return the sample distribution.
#
#    References
#    ----------
#    [1] Marr, C., Geertz, M., Hütt, M.-T., Muskhelishvili, G., 2008.
#        Dissecting the logical types of network control in gene expression profiles.
#        BMC Systems Biology 2, 18.
#    """
#    random_num = int(random_num)
#    if trn is None or trn.size() == 0:
#        return np.nan
#    node2level = {node: lvl for (node, lvl) in izip(active, levels)}
#    (regulators, slaves) = nets.split_regulators(trn)
#    # in TRN structure the out-hubs and spokes differentiation matters
#    reg_levels = [node2level[node] for node in regulators]
#    slave_levels = [node2level[node] for node in slaves]
#    sample = [shuff.continuous_fixed_regulator_sample(trn, regulators, reg_levels,
#            slaves, slave_levels, evaluate=measure) for i in xrange(random_num)]
#    z_score = compute_zscore(measure(trn, node2level), sample)
#    if return_sample:
#        return (z_score, sample)
#    else:
#        return z_score
#
#def delayed_continuous_digital_control(trn, active, levels,
#        delayed_levels, measure=ms.delayed_continuous_functional_coherence):
#    """
#    Compute the continuous digital control in a transcriptional regulatory
#    network (TRN).
#
#    This function computes the continuous digital control using
#    expression levels of nodes in a TRN as a simple null model. Expression
#    levels are considered at two different time points: if there is a link from
#    u to v then the expression level for u at time t is compared with the
#    expression level of v at time point t + 1.
#
#    Parameters
#    ----------
#    trn: TRN or nx.(Multi)DiGraph
#        Effective TRN.
#    active: list
#        Ordered collection of active nodes.
#    levels: list
#        Corresponding expression levels of active nodes (same shape and ordering
#        expected) at time point t.
#    delayed_levels: list
#        Corresponding expression levels of active nodes (same shape and ordering
#        expected) at time point t + 1.
#    measure: callable (optional)
#        Takes the effective network and expression level map and returns the
#        magnitude of control.
#
#    References
#    ----------
#    [1] 
#    """
#    #TODO: complete and test
#    if trn is None or trn.size() == 0:
#        return np.nan
#    node2level = {node: lvl for (node, lvl) in izip(active, levels)}
#    node2delayed = {node: lvl for (node, lvl) in izip(active, delayed_levels)}
#    return measure(trn, node2level, node2delayed)
#
#def delayed_continuous_digital_ctc(trn, active, levels,
#        delayed_levels, random_num=1E04, return_sample=False,
#        measure=ms.delayed_continuous_functional_coherence):
#    """
#    Compute the continuous digital control type confidence (CTC) in a
#    transcriptional regulatory network (TRN).
#
#    This function computes a Z-score of continuous digital control using
#    expression levels of nodes in a TRN as a simple null model. Expression
#    levels are considered at two different time points: if there is a link from
#    u to v then the expression level for u at time t is compared with the
#    expression level of v at time point t + 1.
#
#    Parameters
#    ----------
#    trn: TRN or nx.(Multi)DiGraph
#        Effective TRN.
#    active: list
#        Ordered collection of active nodes.
#    levels: list
#        Corresponding expression levels of active nodes (same shape and ordering
#        expected) at time point t.
#    delayed_levels: list
#        Corresponding expression levels of active nodes (same shape and ordering
#        expected) at time point t + 1.
#    measure: callable (optional)
#        Takes the effective network and expression level map and returns the
#        magnitude of control.
#    random_num: int (optional)
#        Size of the sample distribution of continuous digital control in the
#        null model.
#    return_sample: bool (optional)
#        Whether or not to return the sample distribution.
#
#    References
#    ----------
#    [1] 
#    """
#    #TODO: complete and test
#    random_num = int(random_num)
#    if trn is None or trn.size() == 0:
#        return np.nan
#    node2level = {node: lvl for (node, lvl) in izip(active, levels)}
#    node2delayed = {node: lvl for (node, lvl) in izip(active, delayed_levels)}
#    sample = [shuff.delayed_continuous_sample(trn, active, levels, delayed_levels,
#            evaluate=measure) for i in xrange(random_num)]
#    z_score = compute_zscore(measure(trn, node2level, node2delayed), sample)
#    if return_sample:
#        return (z_score, sample)
#    else:
#        return z_score

def analog_control(effective, measure=ms.discrete_total_ratio):
    """
    Compute the analog control of an effective gene proximity network (GPN).

    Parameters
    ----------
    effective: GPN or nx.(Multi)Graph
        Effective GPN.
    measure: callable (optional)
        Takes the effective network as its only argument and returns the
        magnitude of control within it.

    References
    ----------
    [1] Marr, C., Geertz, M., Hütt, M.-T., Muskhelishvili, G., 2008.
        Dissecting the logical types of network control in gene expression profiles.
        BMC Systems Biology 2, 18.
    """
    if effective is None or effective.size() == 0:
        return np.nan
    return measure(effective)

def analog_ctc(effective, reference, measure=ms.discrete_total_ratio,
        random_num=1E04, return_sample=False):
    """
    Compute the analog control type confidence (CTC) of an effective gene
    proximity network (GPN).

    This function computes a Z-score of analog control using a reference
    (complete) GPN as a simple null model.

    Parameters
    ----------
    effective: GPN or nx.(Multi)Graph
        Effective GPN.
    reference: GPN or nx.(Multi)Graph
        Complete GPN.
    measure: callable (optional)
        Takes the effective network as its only argument and returns the
        magnitude of control within it.
    random_num: int (optional)
        Size of the sample distribution of analog control in the null model.
    return_sample: bool (optional)
        Whether or not to return the sample distribution.

    References
    ----------
    [1] Marr, C., Geertz, M., Hütt, M.-T., Muskhelishvili, G., 2008.
        Dissecting the logical types of network control in gene expression profiles.
        BMC Systems Biology 2, 18.
    """
    random_num = int(random_num)
    if effective is None or effective.size() == 0:
        return np.nan
    size = len(effective)
    sample = [shuff.active_sample(reference, size, evaluate=measure) for i in xrange(random_num)]
    z_score = compute_zscore(measure(effective), sample)
    if return_sample:
        return (z_score, sample)
    else:
        return z_score

#def continuous_analog_control(gpn, active, levels,
#        measure=ms.continuous_abs_coherence):
#    """
#    Compute the continuous analog control of a gene proximity network (GPN).
#
#    Uses expression levels of nodes in a GPN in order to evaluate the magnitude of
#    control.
#
#    Parameters
#    ----------
#    gpn: GPN or nx.(Multi)Graph
#        Effective GPN.
#    active: list
#        Ordered collection of active nodes.
#    levels: list
#        Corresponding expression levels of active nodes (same shape and ordering
#        as active expected).
#    measure: callable (optional)
#        Takes the effective network and expression level map and returns the
#        magnitude of control.
#
#    References
#    ----------
#    [1] 
#    """
#    if gpn is None or gpn.size() == 0:
#        return np.nan
#    node2level = {node: lvl for (node, lvl) in izip(active, levels)}
#    return measure(gpn, node2level)
#
#def continuous_analog_ctc(gpn, active, levels, measure=ms.continuous_abs_coherence,
#        random_num=1E04, return_sample=False):
#    """
#    Compute the continuous analog control type confidence (CTC) in a gene
#    proximity network (GPN).
#
#    This function computes a Z-score of continuous analog control using
#    expression levels of nodes in a GPN as a simple null model.
#
#    Parameters
#    ----------
#    gpn: GPN or nx.(Multi)Graph
#        Effective GPN.
#    active: list
#        Ordered collection of active nodes.
#    levels: list
#        Corresponding expression levels of active nodes (same shape and ordering
#        as active expected).
#    measure: callable (optional)
#        Takes the effective network and expression level map and returns the
#        magnitude of control.
#    random_num: int (optional)
#        Size of the sample distribution of continuous analog control in the
#        null model.
#    return_sample: bool (optional)
#        Whether or not to return the sample distribution.
#
#    References
#    ----------
#    [1] 
#    """
#    random_num = int(random_num)
#    if gpn is None or gpn.size() == 0:
#        return np.nan
#    node2level = {node: lvl for (node, lvl) in izip(active, levels)}
#    sample = [shuff.continuous_sample(gpn, active, levels, evaluate=measure)\
#            for i in xrange(random_num)]
#    z_score = compute_zscore(measure(gpn, node2level), sample)
#    if return_sample:
#        return (z_score, sample)
#    else:
#        return z_score

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
    if rxn_centric is np.nan:
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
        return np.nan
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
    if rxn_centric is np.nan:
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
        return np.nan
    sample = [shuff.active_sample(rxn_centric, size, measure) for i in xrange(int(random_num))]
    z_score = compute_zscore(measure(original), sample)
    if return_sample:
        return (z_score, sample)
    else:
        return z_score

