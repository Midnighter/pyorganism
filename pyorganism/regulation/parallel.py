# -*- coding: utf-8 -*-


"""
===========================================
PyOrganism Regulation Parallel Applications
===========================================

:Authors:
    Moritz Emanuel Beber
:Date:
    2013-02-20
:Copyright:
    Copyright(c) 2013 Jacobs University of Bremen. All rights reserved.
:File:
    parallel.py
"""


__all__ = ["digital_ctc", "continuous_digital_ctc", "analog_ctc",
        "continuous_analog_ctc","metabolic_coherence"]


import logging
import numpy
import random
import re

from itertools import izip
from IPython.parallel import interactive, LoadBalancedView

from . import elements as elem
from .. import miscellaneous as misc
from . import control as cntrl
from . import measures as msr
from ..statistics import compute_zscore
from ..parallel.utils import multi_apply


LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(misc.NullHandler())


def digital_ctc(d_view, trn, active, random_num=1E04, return_sample=False,
        lb_view=None, **kw_args):
    """
    Compute the digital control type confidence of an effective TRN using a list of
    differentially expressed genes.

    Parameters
    ----------
    d_view: `IPython.parallel.DirectView`
        The view is necessary to communicate with the engines.
    trn: nx.(Multi)DiGraph
        Static transcriptional regulatory network.
    active: iterable
        An iterable with gene instances that are differentially expressed.
    lb_view: `IPython.parallel.LoadBalancedView`
        Using a load-balanced view can yield speed improvements.

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
    original = cntrl.setup_trn(trn, active)
    if original is numpy.nan:
        return original
    d_view.execute("import random", block=True)
    d_view.push(dict(network=trn, ratio=msr.discrete_total_ratio), block=True)
    # new null model separates TFs and genes
    t_factors = set(node for node in trn if isinstance(node, elem.TranscriptionFactor))
    genes = set(node for node in trn if isinstance(node, elem.Gene))
    d_view.push(dict(genes=genes, tfs=t_factors), block=True)
    # separate numbers
    tf_num = sum(int(isinstance(item, elem.TranscriptionFactor)) for item in original)
    gene_num = len(original) - tf_num
    d_view.push(dict(gene_num=gene_num, tf_num=tf_num), block=True)
    LOGGER.info("picked %d transcription factors", tf_num)
    if isinstance(lb_view, LoadBalancedView):
        num_krnl = len(lb_view)
        chunk = random_num // num_krnl // 2
        results = multi_apply(lb_view, _trn_sample, random_num, chunksize=chunk)
    else:
        results = multi_apply(d_view, _trn_sample, random_num)
    samples = list(results)
#    LOGGER.info("parallel speed-up was %.3g",
#            results.serial_time / results.wall_time)
    orig_ratio = msr.discrete_total_ratio(original)
    z_score = compute_zscore(orig_ratio, samples)
    if return_sample:
        return (z_score, samples)
    else:
        return z_score

#TODO: adjust to operon based
def continuous_digital_ctc(d_view, trn, active, expr_levels, random_num=1E04,
        return_sample=False, lb_view=None,
        evaluater=msr.continuous_functional_coherence, **kw_args):
    """
    Compute a continuous digital control type confidence of given gene
    expression levels in the effective TRN.

    Parameters
    ----------
    d_view: `IPython.parallel.DirectView`
        The view is necessary to communicate with the engines.
    trn: nx.(Multi)DiGraph
        Static transcriptional regulatory network.
    active: iterable
        An iterable with gene instances.
    expr_levels: iterable
        An iterable in the same order as ``active`` with expression levels.
    lb_view: `IPython.parallel.LoadBalancedView`
        Using a load-balanced view can yield speed improvements.

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
    original = cntrl.setup_trn(trn, active)
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
    (op_net, op2level) = msr._setup_continuous_operon_based(original, orig2level)
    # in TRN structure the out-hubs and spokes differentiation matters
#    out_ops = set(node for (node, deg) in op_net.out_degree_iter() if deg > 0)
#    in_ops = set(op_net.nodes_iter()).difference(out_ops)
#    out_levels = [op2level[op] for op in out_ops]
#    in_levels = [op2level[op] for op in in_ops]
#    d_view.execute("import numpy", block=True)
#    d_view.push(dict(network=op_net, , izip=izip,
#            continuous_coherence=evaluater), block=True) # TODO
#    # TODO pushed active, call remote without parameters
#    if isinstance(lb_view, LoadBalancedView):
#        num_krnl = len(lb_view)
#        chunk = random_num // num_krnl // 2
#        results = lb_view.map(_trn_sample_expression_levels,
#                [active for i in xrange(random_num)],
#                block=False, ordered=False, chunksize=chunk)
#    else:
#        results = d_view.map(_trn_sample_expression_levels,
#                [active for i in xrange(random_num)],
#                block=False)
#    sample = list(results)
#    LOGGER.info("parallel speed-up was %.3g",
#            results.serial_time / results.wall_time)
#    orig_ratio = evaluater(op_net, op2level)
#    z_score = compute_zscore(orig_ratio, sample)
#    if return_sample:
#        return (z_score, sample)
#    else:
#        return z_score

def analog_ctc(d_view, gpn, active, random_num=1E04, return_sample=False,
        lb_view=None, **kw_args):
    """
    Compute the analog control from an effective GPN.

    Parameters
    ----------
    d_view: `IPython.parallel.DirectView`
        The view is necessary to communicate with the engines.
    active: iterable
        An iterable with names of genes that are active in a specific
        condition.
    lb_view: `IPython.parallel.LoadBalancedView`
        Using a load-balanced view can yield speed improvements.

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
    original = cntrl.setup_gpn(gpn, active)
    if original is numpy.nan:
        return original
    d_view.execute("import random", block=True)
    d_view.push(dict(network=gpn, total_ratio=msr.discrete_total_ratio), block=True)
    size = len(original)
    sizes = [size for i in xrange(random_num)]
    if isinstance(lb_view, LoadBalancedView):
        num_krnl = len(lb_view)
        chunk = random_num // num_krnl // 2
        results = lb_view.map(_active_sample, sizes, block=False,
                ordered=False, chunksize=chunk)
    else:
        results = d_view.map(_active_sample, sizes, block=False)
    samples = list(results)
    LOGGER.info("parallel speed-up was %.3g",
            results.serial_time / results.wall_time)
    z_score = compute_zscore(msr.discrete_total_ratio(original), samples)
    if return_sample:
        return (z_score, samples)
    else:
        return z_score

#TODO: adjust to operon based
def continuous_analog_ctc(d_view, organism, active, expr_levels, random_num=1E04,
        return_sample=False, lb_view=None,
        evaluater=msr.continuous_abs_coherence, **kw_args):
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
    original = msr.effective_network(organism.gpn, active)
    size = len(original)
    if size == 0:
        LOGGER.warn("empty effective network")
        return numpy.nan
    d_view.execute("import numpy", block=True)
    rnd_levels = numpy.array(expr_levels, dtype=float, copy=True)
    d_view.push(dict(network=original, expr_levels=rnd_levels, izip=izip,
            continuous_coherence=evaluater), block=True)
    active_genes = [active for i in xrange(random_num)]
    if isinstance(lb_view, LoadBalancedView):
        num_krnl = len(lb_view)
        chunk = random_num // num_krnl // 2
        results = lb_view.map(_gpn_operon_based_sampling, active_genes, block=False,
                ordered=False, chunksize=chunk)
    else:
        results = d_view.map(_gpn_operon_based_sampling, active_genes, block=False)
    sample = list(results)
    gene2level = dict(izip(active, expr_levels))
    orig_ratio = evaluater(original, gene2level)
    z_score = compute_zscore(orig_ratio, sample)
    if return_sample:
        return (z_score, sample)
    else:
        return z_score


def metabolic_coherence(d_view, organism, active, bnumber2gene, rxn_centric=None,
        random_num=1E04, return_sample=False, lb_view=None, **kw_args):
    """
    Compute the metabolic coherence (MC) from an effective metabolic
    network.

    Parameters
    ----------
    d_view: `IPython.parallel.DirectView`
        The view is necessary to communicate with the engines.
    active: iterable
        An iterable with actively expressed genes.
    lb_view: `IPython.parallel.LoadBalancedView`
        Using a load-balanced view can yield speed improvements.

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
        if organism.metabolic_network is None:
            LOGGER.warn("no metabolic network")
            return numpy.nan
        else:
            rxn_centric = organism.metabolic_network.to_reaction_centric()
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
    original = msr.effective_network(rxn_centric, active_reactions)
    size = len(original)
    if size == 0:
        LOGGER.warn("empty effective network")
        return numpy.nan
    d_view.execute("import random", block=True)
    d_view.push(dict(network=rxn_centric, total_ratio=msr.discrete_total_ratio), block=True)
    sizes = [size for i in xrange(int(random_num))]
    if isinstance(lb_view, LoadBalancedView):
        num_krnl = len(lb_view)
        chunk = random_num // num_krnl // 2
        results = lb_view.map(_active_sample, sizes, block=False,
                ordered=False, chunksize=chunk)
    else:
        results = d_view.map(_active_sample, sizes, block=False)
    samples = list(results)
    LOGGER.info("parallel speed-up was %.3g",
            results.serial_time / results.wall_time)
    z_score = compute_zscore(msr.discrete_total_ratio(original), samples)
    if return_sample:
        return (z_score, samples)
    else:
        return z_score

@interactive
def _active_sample(size):
    """
    Sample from affected genes as null model.
    """
    # make use of global variable `network`, `ratio`
    local_net = network
    sample = random.sample(local_net.nodes(), size)
    subnet = local_net.subgraph(sample)
    result = ratio(subnet)
    return result

@interactive
def _trn_sample():
    """
    Sample from affected genes and transcription factors as null model.
    """
    # make use of global variables `network`, `ratio`, `tfs`, `genes`, `tf_num`, `gene_num`
    trn = network
    local_tfs = random.sample(tfs, tf_num)
    local_genes = random.sample(genes, gene_num)
    subnet = trn.subgraph(local_tfs + local_genes)
    result = ratio(subnet)
    return result

#def robustness(self, control_type, active, fraction=0.1,
#        random_num=1E04, control_num=1E04, **kw_args):
#    """
#    Cut a fraction of active genes.
#    """
#    kw_args["random_num"] = control_num
#    size = len(active)
#    cut_num = int(round(size * fraction))
#    LOGGER.info("cutting {0:d}/{1:d} genes".format(cut_num, size))
#    samples = [jackknife(active, cut_num) for i in xrange(int(random_num))]
#    distribution = [control_type(sample, **kw_args) for sample in samples]
#    return distribution
#
#def robustness_with_replacement(self, control_type, active, replacement, fraction=0.1,
#        random_num=1E04, control_num=1E04, **kw_args):
#    """
#    Replace a fraction of the active genes with other genes.
#    """
#    kw_args["random_num"] = control_num
#    size = len(active)
#    replace_num = int(round(size * fraction))
#    LOGGER.info("replacing {0:d}/{1:d} genes".format(replace_num, size))
#    samples = [jack_replace(active, replacement, replace_num)\
#            for i in range(int(random_num))]
#    distribution = [control_type(sample, **kw_args) for sample in samples]
#    return distribution

@interactive
def _trn_sample_expression_levels(active):
    # make use of global variables `network`, `expr_levels`,
    # `continuous_coherence`, `active`
    local_trn = network
#    local_nodes = active
    local_levels = numpy.array(expr_levels, dtype=float, copy=True)
    numpy.random.shuffle(local_levels)
    # build map here because expr_level are shuffled
    gene2level = dict(izip(active, local_levels))
    return continuous_coherence(local_trn, gene2level)

def continuous_trn_operon_sampling(op_net, out_ops, out_levels, in_ops,
        in_levels):
    # select out- and in-ops and then compute similarity based on different
    # null-models
    op2level = dict(izip(out_ops, numpy.random.shuffle(out_levels)))
    op2level.update(izip(in_ops, numpy.random.shuffle(in_levels)))
    return continuous_coherence(op_net, op2level)

@interactive
def _gpn_operon_based_sampling(active):
    # make use of global variables `network`, `expr_levels`,
    # `continuous_coherence`
    local_gpn = network
    local_levels = expr_levels
    all_ops = set(op for gene in active for op in gene.operons)
    orig_gene2level = dict(izip(active, local_levels))
    gene2level = dict()
    no_op = list()
    count = 0
    for gene in active:
        if gene in gene2level:
            continue
        if not gene.operons:
            no_op.append(gene)
            continue
        # multiple operons per gene exist in older versions of RegulonDB
        # we pick shortest of those operons
        ops = list(gene.operons)
        lengths = [len(op.genes) for op in ops]
        op = ops[numpy.argsort(lengths)[0]]
        targets = list(all_ops.difference(set([op])))
        rnd_op = targets[numpy.random.randint(len(targets))]
#        all_ops.difference_update(set([rnd_op]))
        base_level = numpy.nan
        for gn in rnd_op.genes:
            if gn in orig_gene2level:
                base_level = orig_gene2level[gn]
                break
        if numpy.isnan(base_level) and count < 10:
            count += 1
            LOGGER.warn("no change in base level")
        for (i, gn) in enumerate(op.genes):
            gene2level[gn] = base_level
    no_op_levels = [orig_gene2level[gene] for gene in no_op]
    numpy.random.shuffle(no_op_levels)
    gene2level.update(izip(no_op, no_op_levels))
    return continuous_coherence(local_gpn, gene2level)

