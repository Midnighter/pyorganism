# -*- coding: utf-8 -*-


"""
================================
PyOrganism Parallel Applications
================================

:Authors:
    Moritz Emanuel Beber
:Date:
    2013-02-20
:Copyright:
    Copyright(c) 2013 Jacobs University of Bremen. All rights reserved.
:File:
    parallel.py
"""


__all__ = ["digital_ctc", "analog_ctc", "metabolic_coherence"]


import logging
import numpy
import random
import re

from IPython.parallel import interactive, LoadBalancedView

from .regulation import elements as elem
from . import miscellaneous as misc
from .statistics import compute_zscore
from .organism import effective_network, total_ratio


LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(misc.NullHandler())

OPTIONS = misc.OptionsManager.get_instance()


def digital_ctc(d_view, organism, active, random_num=1E04, return_sample=False,
        lb_view=None, **kw_args):
    """
    Compute the digital control type confidence of an effective TRN.

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
    if organism.trn is None or organism.trn.size() == 0:
        LOGGER.warn("empty transcriptional regulatory network")
        return numpy.nan
    active = set([gene.regulatory_product if\
        isinstance(gene.regulatory_product, elem.TranscriptionFactor) else\
        gene for gene in active])
    original = effective_network(organism.trn, active)
    size = len(original)
    if size == 0:
        LOGGER.warn("empty effective network")
        return numpy.nan
    d_view.execute("import random", block=True)
    d_view.push(dict(network=organism.trn, total_ratio=total_ratio), block=True)
    sizes = [size for i in xrange(int(random_num))]
    if isinstance(lb_view, LoadBalancedView):
        num_krnl = len(lb_view)
        chunk = random_num // num_krnl // num_krnl
        results = lb_view.map(_active_sample, sizes, block=False,
                ordered=False, chunksize=chunk)
    else:
        results = d_view.map(_active_sample, sizes, block=False)
    samples = list(results)
    LOGGER.info("parallel speed-up was %.3g",
            results.serial_time / results.wall_time)
    z_score = compute_zscore(total_ratio(original), samples)
    if return_sample:
        return (z_score, samples)
    else:
        return z_score

def analog_ctc(d_view, organism, active, random_num=1E04, return_sample=False,
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
    if organism.gpn is None or organism.gpn.size() == 0:
        LOGGER.warn("empty gene proximity network")
        return numpy.nan
    original = effective_network(organism.gpn, active)
    size = len(original)
    if size == 0:
        LOGGER.warn("empty effective network")
        return numpy.nan
    d_view.execute("import random", block=True)
    d_view.push(dict(network=organism.gpn, total_ratio=total_ratio), block=True)
    sizes = [size for i in xrange(int(random_num))]
    if isinstance(lb_view, LoadBalancedView):
        num_krnl = len(lb_view)
        chunk = random_num // num_krnl // num_krnl
        results = lb_view.map(_active_sample, sizes, block=False,
                ordered=False, chunksize=chunk)
    else:
        results = d_view.map(_active_sample, sizes, block=False)
    samples = list(results)
    LOGGER.info("parallel speed-up was %.3g",
            results.serial_time / results.wall_time)
    z_score = compute_zscore(total_ratio(original), samples)
    if return_sample:
        return (z_score, samples)
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
    original = effective_network(rxn_centric, active_reactions)
    size = len(original)
    if size == 0:
        LOGGER.warn("empty effective network")
        return numpy.nan
    d_view.execute("import random", block=True)
    d_view.push(dict(network=rxn_centric, total_ratio=total_ratio), block=True)
    sizes = [size for i in xrange(int(random_num))]
    if isinstance(lb_view, LoadBalancedView):
        num_krnl = len(lb_view)
        chunk = random_num // num_krnl // num_krnl
        results = lb_view.map(_active_sample, sizes, block=False,
                ordered=False, chunksize=chunk)
    else:
        results = d_view.map(_active_sample, sizes, block=False)
    samples = list(results)
    LOGGER.info("parallel speed-up was %.3g",
            results.serial_time / results.wall_time)
    z_score = compute_zscore(total_ratio(original), samples)
    if return_sample:
        return (z_score, samples)
    else:
        return z_score

@interactive
def _active_sample(size):
    """
    Sample from affected genes as null model.
    """
    # make use of global variable `network`
    local_net = network
    sample = random.sample(local_net.nodes(), size)
    subnet = local_net.subgraph(sample)
    result = total_ratio(subnet)
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

