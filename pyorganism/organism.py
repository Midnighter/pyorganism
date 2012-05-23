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

from operator import itemgetter
from . import miscellaneous as misc
from .parsers import open_file, parser_warning
from .errors import PyOrganismError


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

    def __init__(self, name):
        """
        Parameters
        ----------
        name: str
            The name of the organism. Should be unique but that's not a
            requirement.
        """
        object.__init__(self)
        self.trn = None
        self.gpn = None
        self.go = None
        self.couplons = None
        self.activity = dict()
        self.metabolism = None

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
            if len(partners) > 4 and idn(partners) and partners[2]:
                genes.append((idn(partners), float(partners[2])))
            else:
                warnings(line)
                warnings = LOGGER.warn
        self.activity[description] = genes
        return genes

    def digital_control(self, active):
        """
        Compute the digital control from the effective TRN.

        Parameters
        ----------
        active: iterable
            An iterable with names of genes that are active in a specific
            condition.

        Warning
        -------
        Unknown gene names are silently ignored.
        """
        assert self.trn
        subnet = self.trn.subgraph(active)
        controlled = sum(1 for (node, deg) in subnet.degree_iter() if deg > 0)
        return controlled / float(subnet.order())

    def analog_control(self, active):
        """
        Compute the digital control from the effective GPN.

        Parameters
        ----------
        active: iterable
            An iterable with names of genes that are active in a specific
            condition.

        Warning
        -------
        Unknown gene names are silently ignored.
        """
        assert self.gpn
        subnet = self.gpn.subgraph(active)
        controlled = sum(1 for (node, deg) in subnet.degree_iter() if deg > 0)
        return controlled / float(subnet.order())

