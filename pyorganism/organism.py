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
import multiprocessing

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
        self._pool = multiprocessing.Pool(OPTIONS.num_cpu)

    def __str__(self):
        return str(self.name)

    def __unicode__(self):
        return unicode(self.name)

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

    def digital_control(self, active):
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

    def digital_ctc(self, active, random_num=1E04, parallel=False,
            return_sample=False):
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
        if parallel:
            map_func = self._pool.map
        else:
            map_func = map
        active = set([gene.regulatory_product if type(gene.regulatory_product) ==\
                elem.TranscriptionFactor else gene for gene in active])
        original = effective_network(self.trn, active)
        sample = map_func(active_sample, [(self.trn, len(original))] * int(random_num))
        z_score = compute_zscore(marr_ratio(original), sample)
        if return_sample:
            return (z_score, sample)
        else:
            return z_score

    def analog_control(self, active):
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

    def analog_ctc(self, active, random_num=1E04, parallel=False,
            return_sample=False):
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
        if parallel:
            map_func = self._pool.map
        else:
            map_func = map
        original = effective_network(self.gpn, active)
        sample = map_func(active_sample, [(self.gpn, len(original))] * int(random_num))
        z_score = compute_zscore(marr_ratio(original), sample)
        if return_sample:
            return (z_score, sample)
        else:
            return z_score

    def robustness(self, control_type, active, replacement, fraction=0.1,
            random_num=1E04, parallel=False):
        """
        Replace a fraction of the active genes with other genes.
        """
        if parallel:
            map_func = self._pool.map
        else:
            map_func = map
        kw_args["random_num"] = random_num
        kw_args["parallel"] = parallel
        size = len(active)
        replace_num = int(round(size * fraction))
        LOGGER.info("replacing {0:d}/{1:d} genes".format(replace_num, size))
        samples = map_func(jackknife, [(active, replacement, replace_num)] * int(random_num))
        return distribution


def effective_network(network, active):
    """
    Return the effective network imposed by a subset of nodes.
    """
    subnet = network.subgraph(active)
    size = len(subnet)
    diff = len(active) - size
    LOGGER.info("{0:d} ignored gene(s)".format(diff))
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

def active_sample((network, size)):
    """
    Sample from affected genes as null model.
    """
    sample = random.sample(network.nodes(), size)
    subnet = network.subgraph(sample)
    result = marr_ratio(subnet)
    return result

def jackknife((control_type, active, replacement, replace_num, random_num)):
    positions = random.sample(range(len(active)), replace_num)
    tmp = copy(active)
    replace = random.sample(replacement, replace_num)
    for (i, k) in enumerate(positions):
        tmp[k] = replace[i]
    return tmp

