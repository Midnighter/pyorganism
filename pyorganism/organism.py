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
import warnings
#import multiprocessing

from .regulation import elements as elem
from . import miscellaneous as misc
from .io.generic import open_file, read_tabular
#from .errors import PyOrganismError
from .statistics import compute_zscore


LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(misc.NullHandler())


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
        if self.trn.size() == 0:
            LOGGER.warn("empty transcriptional regulatory network")
        active = set([gene.regulatory_product if type(gene.regulatory_product) ==\
                elem.TranscriptionFactor else gene for gene in active])
        subnet = self.trn.subgraph(active)
        diff = len(active) - len(subnet)
        if diff > 0:
            warnings.warn("{0:d} ignored genes".format(diff))
        if len(subnet) == 0:
            return numpy.nan
        controlled = sum(1 for (node, deg) in subnet.degree_iter() if deg > 0)
#        return controlled / float(subnet.order())
        size = subnet.order()
        if controlled == size:
            return numpy.nan
        else:
            return controlled / float(subnet.order() - controlled)

    def digital_ctc(self, active, random_num=10000, parallel=False,
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
        original = self.digital_control(active)
        length = len(active)
#        if parallel:
#            pool = multiprocessing.Pool()
#            map = pool.map
        sample = map(self.digital_control, [random.sample(self.genes, length)\
                for i in range(random_num)])
        z_score = compute_zscore(original, sample)
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
        if self.gpn.size() == 0:
            LOGGER.warn("empty gene proximity network")
        subnet = self.gpn.subgraph(active)
        diff = len(active) - len(subnet)
        if diff > 0:
            warnings.warn("{0:d} ignored genes".format(diff))
        if len(subnet) == 0:
            return numpy.nan
        controlled = sum(1 for (node, deg) in subnet.degree_iter() if deg > 0)
#        return controlled / float(subnet.order())
        size = subnet.order()
        if controlled == size:
            return numpy.nan
        else:
            return controlled / float(subnet.order() - controlled)

    def analog_ctc(self, active, random_num=10000, parallel=False,
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
        original = self.analog_control(active)
        length = len(active)
#        if parallel:
#            pool = multiprocessing.Pool()
#            map = pool.map
        sample = map(self.analog_control, [random.sample(self.genes, length)\
                for i in range(random_num)])
        z_score = compute_zscore(original, sample)
        if return_sample:
            return (z_score, sample)
        else:
            return z_score

