#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
==============================
PyOrganism Regulatory Networks
==============================

:Authors:
    Moritz Emanuel Beber
:Date:
    2012-05-22
:Copyright:
    Copyright(c) 2012 Jacobs University of Bremen. All rights reserved.
:File:
    networks.py
"""


__all__ = ["TRN", "GPNGenerator"]


import logging
import itertools
import numpy
import networkx as nx

from datetime import date
from operator import itemgetter
from .. import miscellaneous as misc
from ..parsers import open_file, parser_warning
from ..errors import PyOrganismError


LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(misc.NullHandler())

OPTIONS = misc.OptionsManager.get_instance()


class TRN(nx.MultiDiGraph):
    """
    TRN - Transriptional Regulatory Network

    A directed network containing genes and their products (mostly transcription
    factors - TFs). The interactions (lines) must have an attribute denoting the
    regulatory effect:
        * inhibitory: -1
        * activating: 1
        * unknown: 0
        * dual: 2

    Notes
    -----
    Please also read the documentation for networkx.Graph.

    Examples
    --------
    >>> trn = TRN(name="simple")
    >>> trn.add_edge("Ada", "ada", effect=0)
    """

    def read_regulondb(self,
            wsdl="http://regulondb.ccg.unam.mx/webservices/NetWork.jws?wsdl"):
        """
        Retrieve the current version of the TRN from RegulonDB.

        Parameters
        ----------
        wsdl: str (optional)
            The url of the WSDL definition for the DBGET server.

        Notes
        -----
        Requires a SOAPpy installation and an active internet connection.
        """
        # load the required SOAPpy module
        SOAPpy = misc.load_module("SOAPpy", url="http://pywebsvcs.sourceforge.net/")
        # establish connection to DBGET server
        server = SOAPpy.WSDL.Proxy(wsdl)
        interactions = server.getTFGene()
        if not self.name:
            self.name = "TF Gene Network"
        self.graph["date"] = "{:%Y-%m-%d}".format(date.today())
        interactions = interactions.split("\n")
        for line in interactions:
            if line == "":
                continue
            partners = line.split("\t")
            if partners[2] == "repressor":
                self.add_edge(partners[0], partners[1], effect=-1,
                        evidence=partners[3])
            if partners[2] == "activator":
                self.add_edge(partners[0], partners[1], effect=1,
                        evidence=partners[3])
            if partners[2] == "unknown":
                self.add_edge(partners[0], partners[1], effect=0,
                        evidence=partners[3])
            if partners[2] == "dual":
                self.add_edge(partners[0], partners[1], effect=2,
                        evidence=partners[3])

    def read_regulondb_file(self, filename, encoding="utf-8", mode="rb",
            **kw_args):
        """
        Retrieve the TRN from a RegulonDB flat file.

        Parameters
        ----------
        filename: str
            The path to the file containing the interaction information.
        encoding: str (optional)
            The file encoding.
        mode: str (optional)
            The mode used to open a file, should always be read binary.
        """
        if not self.name:
            self.name = filename
        kw_args["encoding"] = encoding
        kw_args["mode"] = mode
        with open_file(filename, **kw_args) as (file_h, ext):
            interactions = file_h.readlines()
        for line in interactions:
            line = line.strip()
            if line == "" or line.startswith("#"):
                continue
            partners = line.split("\t")
            if partners[2] == "-":
                self.add_edge(partners[0], partners[1], effect=-1,
                        evidence=partners[3])
            if partners[2] == "+":
                self.add_edge(partners[0], partners[1], effect=1,
                        evidence=partners[3])
            if partners[2] == "?":
                self.add_edge(partners[0], partners[1], effect=0,
                        evidence=partners[3])
            if partners[2] == "+-":
                self.add_edge(partners[0], partners[1], effect=2,
                        evidence=partners[3])

class GPNGenerator(object):
    """
    GPN - Gene Proximity Network

    An undirected network that contains lines between genes if they are located
    within a certain neighbourhood on the genome.


    Examples
    --------
    """

    def __init__(self):
        """
        Creates a GPNGenerator object.

        A GPNGenerator instance serves the reuse of gene distance information.
        Multiple GPNs can be generated from this information using different
        proximity thresholds.

        Notes
        -----
        `generate_gpn`
        """
        object.__init__(self)
        self.i2name = None
        self.distances = None

    def generate_gpn(self, proximity_threshold, name="", **kw_args):
        """
        Generate the GPN based on previously parsed distance information and a
        given proximity threshold.

        Parameters
        ----------
        proximity_threshold: int
            Maximal distance in base pairs between ORFs on the genome to consider.
        name: str (optional)
            A name for the network.
        Additional keyword arguments are added to the GPN.graph dictionary.

        Notes
        -----
        Please also read the documentation for networkx.Graph.
        """
        proximity_threshold = int(proximity_threshold)
        length = self.distances.shape[0]
        gpn = nx.Graph(name=name, window=proximity_threshold, **kw_args)
        for i in xrange(length - 1):
            for j in xrange(i + 1, length):
                if self.distances[i, j] <= proximity_threshold:
                    gpn.add_edge(self.i2name[i], self.i2name[j])
        return gpn

    def read_regulondb_file(self, filename, gene_identifier="name", comment="#",
            encoding="utf-8", mode="rb", **kw_args):
        """
        Retrieve the gene locations from a RegulonDB flat file and construct the
        GPN using the `proximity_threshold`.

        Parameters
        ----------
        filename: str
            The path to the file containing the interaction information.
        gene_identifier: str (optional)
            Identifiers for the genes used in the network. Can be any of:
                * 'name' for the gene name in that organism
                * 'blattner' for the Blattner number
                * 'regulondb' for the unique identifier assigned by RegulonDB
        comment: str (optional)
            The sign denoting a comment line.
        encoding: str (optional)
            The file encoding.
        mode: str (optional)
            The mode used to open a file, should always be read binary.
        """
        kw_args["encoding"] = encoding
        kw_args["mode"] = mode
        with open_file(filename, **kw_args) as (file_h, ext):
            interactions = file_h.readlines()
        genes = list()
        # choose the column for the gene identifier
        gene_identifier = gene_identifier.lower()
        if gene_identifier == "name":
            idn = itemgetter(1)
        elif gene_identifier == "blattner":
            idn = itemgetter(2)
        elif gene_identifier == "regulondb":
            idn = itemgetter(0)
        else:
            raise PyOrganismError("unrecognised gene identifier '%s'",
                    gene_identifier)
        count = 1
        warnings = parser_warning
        for line in interactions:
            line = line.strip()
            if line == "" or line.startswith(comment):
                continue
            partners = line.split("\t")
            if len(partners) > 4 and idn(partners) and partners[3] and partners[4]:
                name = idn(partners)
                if name == "Phantom Gene":
                    name = "%s %d" % (name, count)
                    count += 1
                # record the identifier, start-, end-position, other information
                genes.append([name, int(partners[3]),
                        int(partners[4])] + partners[5:])
            else:
                warnings(line)
                warnings = LOGGER.warn
        LOGGER.warn("%d phantom genes included", count - 1)
        name = itemgetter(0)
        start = itemgetter(1)
        end = itemgetter(2)
        # assume the maximal end position of the genes is the total length
        genome_length = end(max(genes, key=end))
        length = len(genes)
        self.i2name = dict(itertools.izip(xrange(length),
            (name(gene) for gene in genes)))
        self.distances = numpy.zeros((length, length), dtype=int)
        for i in xrange(length - 1):
            gene_u = genes[i]
            start_u = start(gene_u)
            end_u = end(gene_u)
            for j in xrange(i + 1, length):
                gene_v = genes[j]
                start_v = start(gene_v)
                end_v = end(gene_v)
                # assuming a circular genome here,
                # since RegulonDB is only for E. coli
                if start_u < start_v:
                    diff = start_v - end_u
                else:
                    diff = start_u - end_v
                distance = min(diff, abs(diff - genome_length))
                # we only use the UR triangle of the distances matrix
                self.distances[i, j] = distance

