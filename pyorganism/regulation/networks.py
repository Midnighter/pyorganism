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


__all__ = ["GRN", "TRN", "CouplonGenerator", "GPNGenerator",
    "effective_network", "to_transcription_unit_based", "to_operon_based",
    "trn2grn"]


import logging

import numpy as np
import networkx as nx

from itertools import (product, izip)
from operator import itemgetter, attrgetter
from datetime import date

from .elements import TranscriptionFactor
from ..io.generic import open_file, parser_warning
from .. import miscellaneous as misc
from ..errors import PyOrganismError


LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(misc.NullHandler())

OPTIONS = misc.OptionsManager.get_instance()


class GRN(nx.MultiDiGraph):
    """
    GRN - Gene Regulatory Network

    A directed network containing genes and connections between them if there
    exists a regulatory interaction. Strictly speaking, the interaction is
    between the product of a gene and a target gene but this network holds only
    genes as nodes. Every link must have an attribute denoting the regulatory
    interaction:
        * inhibitory: -1
        * activating: 1
        * unknown: 0
        * dual: 2

    Notes
    -----
    Please also read the documentation for networkx.MultiDiGraph.

    Examples
    --------
    >>> grn = GRN(name="simple")
    >>> grn.add_edge("lacZ", "ada", key=0)
    """

    def __init__(self, data=None, name="", **kw_args):
        super(GRN, self).__init__(data=data, name=name, **kw_args)

    def to_trn(self):
        trn = TRN(name=self.name)
        trn.graph.update(self.graph)
        # do not add nodes since the procedure is ambiguous and
        # trn is defined by regulatory interactions
        for (gene_u, gene_v, k, data) in self.edges_iter(keys=True, data=True):
            if isinstance(gene_u.regulatory_product, TranscriptionFactor):
                u = gene_u.regulatory_product
            else:
                LOGGER.error("gene '%s' has non-transcription factor interaction",
                        gene_u.unique_id)
                continue
            trn.add_edge(u, gene_v, key=k, **data)
            trn.node[u].update(self.node[gene_u])
            trn.node[gene_v].update(self.node[gene_v])
        return trn


class TRN(nx.MultiDiGraph):
    """
    TRN - Transriptional Regulatory Network

    A directed network containing containing transcription factors (TFs) and
    their regulatory targets, i.e., genes. Every link must have an attribute
    denoting the regulatory interaction:
        * inhibitory: -1
        * activating: 1
        * unknown: 0
        * dual: 2

    Notes
    -----
    Please also read the documentation for networkx.MultiDiGraph.

    Examples
    --------
    >>> trn = TRN(name="simple")
    >>> trn.add_edge("Ada", "ada", key=0)
    """

    def __init__(self, data=None, name="", **kw_args):
        super(TRN, self).__init__(data=data, name=name, **kw_args)

    def to_grn(self):
        grn = GRN(name=self.name)
        grn.graph.update(self.graph)
        # do not add nodes since the procedure is ambiguous and
        # trn is defined by regulatory interactions
        for (tf, gene, k, data) in self.edges_iter(keys=True, data=True):
            for src_gene in tf.coded_from:
                grn.add_edge(src_gene, gene, key=k, **data)
                grn.node[src_gene].update(self.node[tf])
                grn.node[gene].update(self.node[gene])
        return grn

    def to_couplons(self, sf_links):
        couplon_gen = CouplonGenerator(self)
        couplon_gen.add_edges_from(sf_links)
        return couplon_gen


class CouplonGenerator(nx.DiGraph):
    """
    """
    def __init__(self, trn, name="", **kw_args):
        """
        Creates a CouplonGenerator instance.

        The CouplonGenerator should be initialised with a TRN. `trn` is a
        keyword argument here to allow for the networkx subgraphing utility.

        Notes
        -----
        `generate_couplon`
        """
        super(CouplonGenerator, self).__init__(data=trn, name=name, **kw_args)

    def generate_couplon(self, nap, sigma_factor):
        """
        Using a network that includes sigma factor, transcription factor, and
        gene interactions, couplons are generated using a nucleoid associated
        protein (NAP) and sigma factor pair to generate subgraphs from the pair
        and their common successors.

        Parameters
        ----------
        nap: str
            A NAP name which can be found among the TFs.
        sigma_factor: str
            A sigma factor.
        """
        # find NAP successors
        nap_targets = set(self.successors_iter(nap))
        # find sigma factor successors
        sigma_targets = set(self.successors_iter(sigma_factor))
        # find the common regulated targets
        common = nap_targets.intersection(sigma_targets)
        # add NAP and sigma factor
        common.add(nap)
        common.add(sigma_factor)
        # return couplon
        couplon = self.subgraph(common)
        couplon.name = "%s-%s Couplon" % (nap, sigma_factor)
        return couplon

    def lookup_regulondb(self, name2gene=None, sep="\t",
            wsdl="http://regulondb.ccg.unam.mx/webservices/NetWork.jws?wsdl"):
        """
        Retrieve the current version of the sigma factor and gene interactions
        from RegulonDB.

        Parameters
        ----------
        name2gene: dict (optional)
            The genes in the interactions between sigma factors and genes parsed
            from RegulonDB are given by their name. If the genes in the TRN are
            identified by their Blattner number of RegulonDB identifier, a
            mapping is required.
        sep: str (optional)
            The column separator; not likely to change.
        wsdl: str (optional)
            The url of the WSDL definition for the DBGET server.

        Notes
        -----
        Requires a SOAPpy installation and an active internet connection.
        """
        def use_mapping(tf, gene, effct, evdnc):
            self.add_edge(tf, name2gene[gene], interaction=effct, evidence=evdnc)

        def no_mapping(tf, gene, effct, evdnc):
            self.add_edge(tf, gene, interaction=effct, evidence=evdnc)

        # load the required SOAPpy module
        SOAPpy = misc.load_module("SOAPpy", url="http://pywebsvcs.sourceforge.net/")
        # establish connection to DBGET server
        server = SOAPpy.WSDL.Proxy(wsdl)
        interactions = server.getTFGene()
        interactions = interactions.split("\n")
        self.graph["date"] = "{:%Y-%m-%d}".format(date.today())
        if not self.name:
            self.name = "Sigma Factor - Gene Network"
        if name2gene:
            add_link = use_mapping
        else:
            add_link = no_mapping
        warnings = parser_warning
        # the last line is an empty one
        for line in interactions[:-1]:
            line = line.strip()
            partners = line.split(sep)
            # sigma factors only have activating roles
            if partners[2] == "activator":
                add_link(partners[0], partners[1], 1, partners[3])
            else:
                warnings(line)
                warnings = LOGGER.warn

    def read_regulondb(self, filename, name2gene=None, sep="\t",
            encoding="utf-8", mode="rb", **kw_args):
        """
        Retrieve sigma factor and gene interactions from a RegulonDB flat file.

        Parameters
        ----------
        filename: str
            The path to the file containing the interaction information.
        name2gene: dict (optional)
            The genes in the interactions between sigma factors and genes parsed
            from RegulonDB are given by their name. If the genes in the TRN are
            identified by their Blattner number of RegulonDB identifier, a
            mapping is required.
        sep: str (optional)
            The column separator; not likely to change.
        encoding: str (optional)
            The file encoding.
        mode: str (optional)
            The mode used to open a file, should always be read binary.
        """
        def use_mapping(tf, gene, effct, evdnc):
            self.add_edge(tf, name2gene[gene], interaction=effct, evidence=evdnc)

        def no_mapping(tf, gene, effct, evdnc):
            self.add_edge(tf, gene, interaction=effct, evidence=evdnc)

        if not self.name:
            self.name = filename
        kw_args["encoding"] = encoding
        kw_args["mode"] = mode
        with open_file(filename, **kw_args) as (file_h, ext):
            interactions = file_h.readlines()
        if name2gene:
            add_link = use_mapping
        else:
            add_link = no_mapping
        warnings = parser_warning
        for line in interactions:
            line = line.strip()
            if line == "" or line.startswith("#"):
                continue
            partners = line.split(sep)
            # sigma factors only have activating roles
            if partners[2] == "+":
                add_link(partners[0], partners[1], 1, partners[3])
            else:
                warnings(line)
                warnings = LOGGER.warn


class GPNGenerator(object):
    """
    GPN - Gene Proximity Network

    An undirected network that contains lines between genes if they are located
    within a certain neighbourhood on the genome.


    Examples
    --------
    """

    def __init__(self, **kw_args):
        """
        Creates a GPNGenerator object.

        A GPNGenerator instance serves the reuse of gene distance information.
        Multiple GPNs can be generated from this information using different
        proximity thresholds.

        Notes
        -----
        `generate_gpn`
        """
        super(GPNGenerator, self).__init__(**kw_args)
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
        gpn.add_nodes_from(self.i2name.itervalues())
        valid = self.distances > -1
        for i in xrange(length - 1):
            for j in xrange(i + 1, length):
                if valid[i, j] and self.distances[i, j] <= proximity_threshold:
                    gpn.add_edge(self.i2name[i], self.i2name[j])
        return gpn

    def parse_genes(self, genes):
        end = attrgetter("position_end")
        # assume the maximal end position of the genes is the total length
        genome_length = end(max(genes, key=end))
        length = len(genes)
        self.i2name = dict(izip(xrange(length), genes))
        self.distances = np.zeros((length, length), dtype=int)
        self.distances.fill(-1)
        diffs = np.zeros(4, dtype=int)
        no_position = set()
        for i in xrange(length - 1):
            gene_u = genes[i]
            if gene_u.position is None or None in gene_u.position:
                no_position.add(gene_u)
                continue
            for j in xrange(i + 1, length):
                gene_v = genes[j]
                if gene_v.position is None or None in gene_v.position:
                    no_position.add(gene_v)
                    continue
                # assuming a circular genome here,
                # since RegulonDB is only for E. coli
                # compute difference between start and end points (with overlap)
                for (k, pair) in enumerate(product(gene_u.position,
                        gene_v.position)):
                    diff = abs(pair[0] - pair[1])
                    diffs[k] = min(diff, genome_length - diff)
                # we only use the UR triangle of the distances matrix
                self.distances[i, j] = diffs.min()
        for gene in no_position:
            LOGGER.warn("no position information for gene '{0}'".format(gene))

    def read_regulondb(self, filename, gene_identifier="name", sep="\t",
            comment="#", encoding="utf-8", mode="rb", **kw_args):
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
        sep: str (optional)
            The column separator; not likely to change.
        comment: str (optional)
            The sign denoting a comment line.
        encoding: str (optional)
            The file encoding.
        mode: str (optional)
            The mode used to open a file, should always be read binary.
        """
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
        # read information from the file
        kw_args["encoding"] = encoding
        kw_args["mode"] = mode
        with open_file(filename, **kw_args) as (file_h, ext):
            information = file_h.readlines()
        # parse the information
        genes = list()
        phantom_num = 0
        warnings = parser_warning
        for line in information:
            line = line.strip()
            if line == "" or line.startswith(comment):
                continue
            partners = line.split(sep)
            if len(partners) > 4 and idn(partners) and partners[3] and partners[4]:
                name = idn(partners)
                if name == "Phantom Gene":
                    phantom_num += 1
                    name = "%s %d" % (name, phantom_num)
                # record the identifier, start-, end-position, other information
                genes.append([name, int(partners[3]),
                        int(partners[4])] + partners[5:])
            else:
                warnings(line)
                warnings = LOGGER.warn
        LOGGER.warn("%d phantom genes included", phantom_num)
        name = itemgetter(0)
        start = itemgetter(1)
        end = itemgetter(2)
        # assume the maximal end position of the genes is the total length
        genome_length = end(max(genes, key=end))
        length = len(genes)
        self.i2name = dict(izip(xrange(length),
            (name(gene) for gene in genes)))
        self.distances = np.zeros((length, length), dtype=int)
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
                # compute difference between start and end points (overlap)
                diff_1 = abs(start_u - end_v)
                diff_1 = min(diff_1, genome_length - diff_1)
                diff_2 = abs(start_v - end_u)
                diff_2 = min(diff_2, genome_length - diff_2)
                diff_3 = abs(start_u - start_v)
                diff_3 = min(diff_3, genome_length - diff_3)
                diff_4 = abs(end_v - end_u)
                diff_4 = min(diff_4, genome_length - diff_4)
                # we only use the UR triangle of the distances matrix
                self.distances[i, j] = min(diff_1, diff_2, diff_3, diff_4)


def effective_network(network, active):
    """
    Return the effective network imposed by a subset of nodes.

    Warning
    -------
    Node instances not in the network are silently ignored.

    """
    subnet = network.subgraph(active)
    size = len(subnet)
    LOGGER.info("{0:d}/{1:d} node(s) effective network - {2:d} entities ignored"\
            .format(size, len(network), len(active) - size))
    return subnet

def split_regulators(net):
    """
    Return a pair of all controlling nodes and others.
    """
    regulators = {node for (node, deg) in net.out_degree_iter() if deg > 0}
    if len(regulators) == 0:
        LOGGER.error("no regulating hubs in TRN")
    slaves = set(net.nodes_iter()) - regulators
    if len(slaves) == 0:
        LOGGER.error("no regulated targets in TRN")
    return (list(regulators), list(slaves))

def to_transcription_unit_based(network):
    tu_net = type(network)()
    for node in network:
        ops = node.get_transcription_units()
        if len(ops) == 0:
            tu_net.add_node(node)
        else:
            tu_net.add_nodes_from(ops)
    if tu_net.is_multigraph():
        for (u, v, inter) in network.edges_iter(keys=True):
            u_ops = u.get_transcription_units()
            if len(u_ops) == 0:
                u_ops = [u]
            v_ops = v.get_transcription_units()
            if len(v_ops) == 0:
                v_ops = [v]
            for (op_1, op_2) in product(u_ops, v_ops):
                tu_net.add_edge(op_1, op_2, key=inter)
    else:
        for (u, v) in network.edges_iter():
            u_ops = u.get_transcription_units()
            if len(u_ops) == 0:
                u_ops = [u]
            v_ops = v.get_transcription_units()
            if len(v_ops) == 0:
                v_ops = [v]
            for (op_1, op_2) in product(u_ops, v_ops):
                tu_net.add_edge(op_1, op_2)
    return tu_net

def to_operon_based(network):
    op_net = type(network)()
    for node in network:
        ops = node.get_operons()
        if len(ops) == 0:
            op_net.add_node(node)
        else:
            op_net.add_nodes_from(ops)
    if op_net.is_multigraph():
        for (u, v, inter) in network.edges_iter(keys=True):
            u_ops = u.get_operons()
            if len(u_ops) == 0:
                u_ops = [u]
            v_ops = v.get_operons()
            if len(v_ops) == 0:
                v_ops = [v]
            for (op_1, op_2) in product(u_ops, v_ops):
                op_net.add_edge(op_1, op_2, key=inter)
    else:
        for (u, v) in network.edges_iter():
            u_ops = u.get_operons()
            if len(u_ops) == 0:
                u_ops = [u]
            v_ops = v.get_operons()
            if len(v_ops) == 0:
                v_ops = [v]
            for (op_1, op_2) in product(u_ops, v_ops):
                op_net.add_edge(op_1, op_2)
    return op_net

def trn2grn(network):
    grn = type(network)()
    grn.graph.update(network.graph)
    # do not add nodes since the procedure is ambiguous and
    # trn is defined by regulatory interactions
    if network.is_multigraph():
        for (tf, gene, k, data) in network.edges_iter(keys=True, data=True):
            for src_gene in tf.coded_from:
                grn.add_edge(src_gene, gene, key=k, **data)
                grn.node[src_gene].update(network.node[tf])
                grn.node[gene].update(network.node[gene])
    else:
        for (tf, gene, data) in network.edges_iter(data=True):
            for src_gene in tf.coded_from:
                grn.add_edge(src_gene, gene, **data)
                grn.node[src_gene].update(network.node[tf])
                grn.node[gene].update(network.node[gene])
    return grn

def setup_metabolic(metabolic, rxn_centric):
    if rxn_centric is None:
        if metabolic is None:
            LOGGER.error("no metabolic network")
            return np.nan
        else:
            rxn_centric = metabolic.to_reaction_centric()
    if rxn_centric.size() == 0:
        LOGGER.error("empty metabolic network")
        return np.nan

