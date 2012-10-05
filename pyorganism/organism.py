#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
==========
PyOrganism
==========

:Authors:
    Moritz Emanuel Beber
:Date:
    2012-05-22
:Copyright:
    Copyright(c) 2012 Jacobs University of Bremen. All rights reserved.
:File:
    organism.py
"""


__all__ = ["Organism"]


import logging
import numpy
import random
import re

from operator import attrgetter
from copy import copy
from .regulation import elements as elem
from . import miscellaneous as misc
from .io.generic import open_file, read_tabular
#from .errors import PyOrganismError
from .statistics import compute_zscore


LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(misc.NullHandler())

OPTIONS = misc.OptionsManager.get_instance()


class Organism(object):
    """
    A representation of a living organism with multiple layers of organisation.

    As many layers of organisation as are available or desired may be included
    in the `Organism` object.

    Notes
    -----

    Examples
    --------


    """

    def __init__(self, name, **kw_args):
        """
        Parameters
        ----------
        name: str
            The name of the organism. Should be unique but that's not a
            requirement.
        """
        super(Organism, self).__init__(**kw_args)
        self.name = name
        self.genes = None
        self.trn = None
        self.gpn = None
        self.go = None
        self.couplons = None
        self.activity = dict()
        self.metabolism = None
        self.metabolic_network = None

    def __str__(self):
        return str(self.name)

    def __unicode__(self):
        return unicode(self.name)

    def read_expression_ration(self, filename, description, sep="\t",
            skip_header=0, comment="#",
            encoding="utf-8", mode="rb", **kw_args):
        """
        Retrieve gene activity from prepared micro array data.

        Parameters
        ----------
        filename: str
            The path to the file containing the interaction information.
        description: str
            A unique identifier for the conditions of the gene activity pattern
            to parse.
        sep: str (optional)
            The separator distinguishing columns in the file.
        gene_identifier: str (optional)
            Identifiers for the genes used in the network. Can be any of:
                * 'name' for the gene name in that organism
                * 'blattner' for the Blattner number
        skip_header: int (optional)
            The number of header lines to skip.
        comment: str (optional)
            The sign denoting a comment line.
        encoding: str (optional)
            The file encoding.
        mode: str (optional)
            The mode used to open a file, should always be read binary.
        """
        kw_args["encoding"] = encoding
        kw_args["mode"] = mode
        genes = list()
        with open_file(filename, **kw_args) as (file_h, ext):
            lines = file_h.readlines()
            for row in read_tabular(lines[skip_header:], sep=sep, comment=comment):
                if len(row) > 1 and row[1]:
                    if row[0]:
                        genes.append((row[0], float(row[1])))
                    else:
                        LOGGER.warn(row)
        self.activity[description] = genes
        return genes

    def read_activity(self, filename, description, sep="\t",
            skip_header=1, comment="#",
            encoding="utf-8", mode="rb", **kw_args):
        """
        Retrieve gene activity from prepared micro array data.

        Parameters
        ----------
        filename: str
            The path to the file containing the interaction information.
        description: str
            A unique identifier for the conditions of the gene activity pattern
            to parse.
        sep: str (optional)
            The separator distinguishing columns in the file.
        gene_identifier: str (optional)
            Identifiers for the genes used in the network. Can be any of:
                * 'name' for the gene name in that organism
                * 'blattner' for the Blattner number
        skip_header: int (optional)
            The number of header lines to skip.
        comment: str (optional)
            The sign denoting a comment line.
        encoding: str (optional)
            The file encoding.
        mode: str (optional)
            The mode used to open a file, should always be read binary.
        """
        kw_args["encoding"] = encoding
        kw_args["mode"] = mode
        genes = list()
        with open_file(filename, **kw_args) as (file_h, ext):
            lines = file_h.readlines()
            for row in read_tabular(lines[skip_header:], sep=sep, comment=comment):
                if len(row) > 3 and row[2]:
                    if row[0]:
                        genes.append((row[0], float(row[2])))
                    elif row[1]:
                        genes.append((row[1], float(row[2])))
                    else:
                        LOGGER.warn(row)
        self.activity[description] = genes
        return genes

    def digital_control(self, active, **kw_args):
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
        active = set([gene.regulatory_product if type(gene.regulatory_product) ==\
                elem.TranscriptionFactor else gene for gene in active])
        subnet = effective_network(self.trn, active)
        return marr_ratio(subnet)

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
        active = set([gene.regulatory_product if type(gene.regulatory_product) ==\
                elem.TranscriptionFactor else gene for gene in active])
        original = effective_network(self.trn, active)
        size = len(original)
        sample = [active_sample(self.trn, size) for i in range(int(random_num))]
        z_score = compute_zscore(marr_ratio(original), sample)
        if return_sample:
            return (z_score, sample)
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
        return marr_ratio(subnet)

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
        sample = [active_sample(self.gpn, size) for i in range(int(random_num))]
        z_score = compute_zscore(marr_ratio(original), sample)
        if return_sample:
            return (z_score, sample)
        else:
            return z_score

    def metabolic_coherence_ratio(self, active, bnumber2gene, **kw_args):
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
        if self.metabolic_network is None or self.metabolic_network.size() == 0:
            LOGGER.warn("empty metabolic network")
            return numpy.nan
        rxn_centric = self.metabolic_network.to_reaction_centric()
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
        return total_ratio(subnet)

    def metabolic_coherence(self, active, bnumber2gene, random_num=1E04,
            return_sample=False, **kw_args):
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
        if self.metabolic_network is None or self.metabolic_network.size() == 0:
            LOGGER.warn("empty metabolic network")
            return numpy.nan
        rxn_centric = self.metabolic_network.to_reaction_centric()
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
        sample = [active_sample(rxn_centric, size) for i in range(int(random_num))]
        z_score = compute_zscore(marr_ratio(original), sample)
        if return_sample:
            return (z_score, sample)
        else:
            return z_score

    def robustness(self, control_type, active, replacement, fraction=0.1,
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
    diff = len(active) - size
    LOGGER.info("{0:d} ignored node(s)".format(diff))
    if size == 0:
        LOGGER.warn("empty effective network")
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
    result = marr_ratio(subnet)
    return result

def jack_replace(active, replacement, replace_num):
    positions = random.sample(range(len(active)), replace_num)
    tmp = copy(active)
    replace = random.sample(replacement, replace_num)
    for (i, k) in enumerate(positions):
        tmp[k] = replace[i]
    return tmp

def jackknife(active, remove_num):
    positions = random.sample(range(len(active)), remove_num)
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

