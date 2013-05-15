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


__all__ = ["digital_control", "digital_ctc", "analog_control", "analog_ctc",
        "continuous_analog_ctc", "metabolic_coherence_ratio", "metabolic_coherence"]


import logging
import numpy
import random
import re

from operator import attrgetter
from copy import copy
from itertools import izip

from . import elements as elem
from .. import miscellaneous as misc
#from ..errors import PyOrganismError
from ..statistics import compute_zscore


LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(misc.NullHandler())

OPTIONS = misc.OptionsManager.get_instance()


def digital_control(organism, active, **kw_args):
    """
    Compute the digital control from an effective TRN.

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
    if self.trn is None or self.trn.size() == 0:
        LOGGER.warn("empty transcriptional regulatory network")
        return numpy.nan
    active = set([gene.regulatory_product if\
        isinstance(gene.regulatory_product, elem.TranscriptionFactor) else\
        gene for gene in active])
    subnet = effective_network(self.trn, active)
    if len(subnet) == 0:
        LOGGER.warn("empty effective network")
        return numpy.nan
    return total_ratio(subnet)

def digital_ctc(self, active, random_num=1E04, return_sample=False, **kw_args):
    """
    Compute the digital control type confidence of an effective TRN.

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
    if self.trn is None or self.trn.size() == 0:
        LOGGER.warn("empty transcriptional regulatory network")
        return numpy.nan
    active = set([gene.regulatory_product if\
        isinstance(gene.regulatory_product, elem.TranscriptionFactor) else\
        gene for gene in active])
    original = effective_network(self.trn, active)
    size = len(original)
    if size == 0:
        LOGGER.warn("empty effective network")
        return numpy.nan
    sample = [active_sample(self.trn, size) for i in range(int(random_num))]
    z_score = compute_zscore(total_ratio(original), sample)
    if return_sample:
        return (z_score, sample)
    else:
        return z_score

def digital_ctc2(self, active, random_num=1E04, return_sample=False, **kw_args):
    random_num = int(random_num)
    if self.trn is None or self.trn.size() == 0:
        LOGGER.warn("empty transcriptional regulatory network")
        return numpy.nan
    active = set([gene.regulatory_product if\
        isinstance(gene.regulatory_product, elem.TranscriptionFactor) else\
        gene for gene in active])
    LOGGER.info("picked %d transcription factors", sum(1 for item in active if
        isinstance(item, elem.TranscriptionFactor)))
    original = effective_network(self.trn, active)
    size = len(original)
    if size == 0:
        LOGGER.warn("empty effective network")
        return numpy.nan
    # new null model separates TFs and genes
    t_factors = set(node for node in self.trn if isinstance(node, elem.TranscriptionFactor))
    genes = set(node for node in self.trn if isinstance(node, elem.Gene))
    # separate numbers
    tf_num = sum(1 for item in original if isinstance(item, elem.TranscriptionFactor))
    gene_num = sum(1 for item in original if isinstance(item, elem.Gene))
    samples = [trn_sample(self.trn, t_factors, tf_num, genes, gene_num) for i in range(int(random_num))]
    z_score = compute_zscore(total_ratio(original), samples)
    if return_sample:
        return (z_score, samples)
    else:
        return z_score

def analog_control(self, active, **kw_args):
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
    if self.gpn is None or self.gpn.size() == 0:
        LOGGER.warn("empty gene proximity network")
        return numpy.nan
    subnet = effective_network(self.gpn, active)
    if len(subnet) == 0:
        LOGGER.warn("empty effective network")
        return numpy.nan
    return total_ratio(subnet)

def analog_ctc(self, active, random_num=1E04, return_sample=False, **kw_args):
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
    if self.gpn is None or self.gpn.size() == 0:
        LOGGER.warn("empty gene proximity network")
        return numpy.nan
    original = effective_network(self.gpn, active)
    size = len(original)
    if size == 0:
        LOGGER.warn("empty effective network")
        return numpy.nan
    sample = [active_sample(self.gpn, size) for i in range(int(random_num))]
    z_score = compute_zscore(total_ratio(original), sample)
    if return_sample:
        return (z_score, sample)
    else:
        return z_score

def continuous_analog_ctc(organism, active, expr_levels, random_num=1E04,
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
    if organism.gpn is None or organism.gpn.size() == 0:
        LOGGER.warn("empty gene proximity network")
        return numpy.nan
    original = effective_network(organism.gpn, active)
    size = len(original)
    if size == 0:
        LOGGER.warn("empty effective network")
        return numpy.nan
    rnd_levels = numpy.array(expr_levels, dtype=float, copy=True)
    sample = [gpn_sample_expression_levels(original, active, rnd_levels)\
            for i in range(random_num)]
    gene2level = dict(izip(active, expr_levels))
    orig_ratio = gpn_expression_level_similarity(original, gene2level)
    z_score = compute_zscore(orig_ratio, sample)
    if return_sample:
        return (z_score, sample)
    else:
        return z_score

def metabolic_coherence_ratio(self, active, bnumber2gene, rxn_centric=None,
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
    subnet = effective_network(rxn_centric, active_reactions)
    if len(subnet) == 0:
        LOGGER.warn("empty effective network")
        return numpy.nan
    return total_ratio(subnet)

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
    z_score = compute_zscore(total_ratio(original), sample)
    if return_sample:
        return (z_score, sample)
    else:
        return z_score

def robustness(self, control_type, active, fraction=0.1,
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

def robustness_with_replacement(self, control_type, active, replacement, fraction=0.1,
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

def stat_expression_forks(self, statistic, active, ori=(3923657, 3924034),
        ter=None):
    end = attrgetter("position_end")
    # assume the maximal end position of the genes is the total length
    genome_length = end(max(self.genes, key=end))
    if ter is None:
        LOGGER.debug("genome length = {0:d}".format(genome_length))
        ter = int(round(genome_length / 2.0 + (ori[1] + ori[0]) / 2.0)) % genome_length
        LOGGER.debug("terminus location = {0:d}".format(ter))
        ter = (ter, ter)
    dist = abs(ori[1] - ter[0])
    full = genome_length - ori[1]
    fork = list()
    counter_fork = list()
    for gene in active:
        if None in gene.position:
            continue
        if gene.position_start > ori[1] and gene.position_start < (ori[1] + dist):
            location = gene.position_start - ori[1]
            fork.append((location, statistic(gene, True)))
        elif gene.position_start < ter[0] and gene.position_start > (ter[0] - dist):
            location = full + gene.position_start
            fork.append((location, statistic(gene, True)))
        else:
            location = ori[0] - gene.position_end
            counter_fork.append((location, statistic(gene, False)))
    return ([x for x in fork if not x[1] is None],
            [x for x in counter_fork if not x[1] is None])


def effective_network(network, active):
    """
    Return the effective network imposed by a subset of nodes.
    """
    subnet = network.subgraph(active)
    size = len(subnet)
    LOGGER.info("{0:d}/{1:d} node(s) effective network - {2:d} entities ignored"\
            .format(size, len(network), len(active) - size))
    return subnet

def marr_ratio(network):
    control = sum(1 for (node, deg) in network.degree_iter() if deg > 0)
    if control == len(network):
        return numpy.nan
    else:
        return control / float(len(network) - control)

def total_ratio(network):
    control = sum(1 for (node, deg) in network.degree_iter() if deg > 0)
    return control / float(len(network))

def active_sample(network, size):
    """
    Sample from affected genes as null model.
    """
    sample = random.sample(network.nodes(), size)
    subnet = network.subgraph(sample)
    result = total_ratio(subnet)
    return result

def trn_sample(trn, tfs, tf_num, genes, gene_num):
    """
    Sample from affected genes and transcription factors as null model.
    """
    local_tfs = random.sample(tfs, tf_num)
    local_genes = random.sample(genes, gene_num)
    subnet = trn.subgraph(local_tfs + local_genes)
    result = total_ratio(subnet)
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

def leading_gc(gene, clockwise):
    if clockwise:
        if gene.strand == "forward":
            return gene.gc_content
    else:
        if gene.strand == "reverse":
            return gene.gc_content

def lagging_gc(gene, clockwise):
    if clockwise:
        if gene.strand == "reverse":
            return gene.gc_content
    else:
        if gene.strand == "forward":
            return gene.gc_content

def gpn_sample_expression_levels(network, active, expr_level):
    numpy.random.shuffle(expr_level)
    # build map here because expr_level are shuffled
    gene2level = dict(izip(active, expr_level))
    return gpn_expression_level_similarity(network, gene2level)

def gpn_expression_level_similarity(network, gene2level):
    return sum(1.0 - abs(gene2level[u] - gene2level[v])\
            for (u, v) in network.edges_iter())
