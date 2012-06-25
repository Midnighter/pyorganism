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

from operator import itemgetter
from . import miscellaneous as misc
from .io.generic import open_file, parser_warning
from .errors import PyOrganismError
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

    def __init__(self, name, *args, **kw_args):
        """
        Parameters
        ----------
        name: str
            The name of the organism. Should be unique but that's not a
            requirement.
        """
        super(Organism, self).__init__(*args, **kw_args)
        self.name = name
        self.trn = None
        self.gpn = None
        self.go = None
        self.couplons = None
        self.activity = dict()
        self.metabolism = None

    def __str__(self):
        return str(self)

    def __unicode__(self):
        return unicode(self)

    def read_activity(self, filename, description, sep="\t",
            gene_identifier="name", skip_header=1, comment="#",
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
        # choose the column for the gene identifier
        gene_identifier = gene_identifier.lower()
        if gene_identifier == "name":
            idn = itemgetter(0)
        elif gene_identifier == "blattner":
            idn = itemgetter(1)
        else:
            raise PyOrganismError("unrecognised gene identifier '%s'",
                    gene_identifier)
        with open_file(filename, **kw_args) as (file_h, ext):
            interactions = file_h.readlines()
        genes = list()
        warnings = parser_warning
        for line in interactions[skip_header:]:
            line = line.strip()
            if line == "" or line.startswith(comment):
                continue
            partners = line.split(sep)
            # functional description is not important here
            if len(partners) > 3 and idn(partners) and partners[2]:
                genes.append((idn(partners), float(partners[2])))
            else:
                warnings(line)
                warnings = LOGGER.warn
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
        subnet = self.trn.subgraph(active)
        if subnet.order() == 0:
            LOGGER.warn("set of active genes unknown")
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
        genes = self.trn.nodes()
        length = len(active)
#        if parallel:
#            pool = multiprocessing.Pool()
#            map = pool.map
        sample = map(self.digital_control, [random.sample(genes, length)\
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
        subnet = self.gpn.subgraph(active)
        if subnet.order() == 0:
            LOGGER.warn("set of active genes unknown")
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
        genes = self.gpn.nodes()
        length = len(active)
#        if parallel:
#            pool = multiprocessing.Pool()
#            map = pool.map
        sample = map(self.analog_control, [random.sample(genes, length)\
                for i in range(random_num)])
        z_score = compute_zscore(original, sample)
        if return_sample:
            return (z_score, sample)
        else:
            return z_score

