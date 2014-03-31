#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
=================================
Metabolic Network Representations
=================================

:Authors:
    Moritz Emanuel Beber
:Date:
    2011-04-13
:Copyright:
    Copyright(c) 2011 Jacobs University of Bremen. All rights reserved.
:File:
    networks.py
"""


__all__ = ["MetabolicNetwork",
        "CompoundCentricNetwork", "CompoundCentricMultiNetwork",
        "ReactionCentricNetwork", "ReactionCentricMultiNetwork",
        "read_edgelist", "write_edgelist"]

import logging
import itertools
import networkx as nx

from collections import defaultdict
from ..errors import PyOrganismError
from .. import miscellaneous as misc
from ..io.generic import open_file, read_tabular
from . import elements as pymet


LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(misc.NullHandler())

OPTIONS = misc.OptionsManager.get_instance()


class MetabolicNetwork(nx.DiGraph):
    """
    """

    def __init__(self, data=None, name="", compartments=None, **kw_args):
        """
        """
        super(MetabolicNetwork, self).__init__(data=data, name=name, **kw_args)
        self.compartments = misc.convert(compartments, set, set())
        self.compounds = set()
        self.reactions = set()
        self.rpairs = None

    def add_edge(self, u, v, **kw_args):
        """
        """
        if isinstance(u, pymet.BasicReaction):
            self.reactions.add(u)
        elif isinstance(u, pymet.BasicCompound):
            self.compounds.add(u)
        else:
            raise TypeError("unidentified metabolic type '{0}'".format(type(u)))
        if isinstance(v, pymet.BasicReaction):
            self.reactions.add(v)
        elif isinstance(v, pymet.BasicCompound):
            self.compounds.add(v)
        else:
            raise TypeError("unidentified metabolic type '{0}'".format(type(v)))
        super(MetabolicNetwork, self).add_edge(u, v, **kw_args)

    def remove_edge(self, u, v):
        """
        """
        super(MetabolicNetwork, self).remove_edge(u, v)
        # removing a link does not necessitate removing the involved nodes
#        if u not in self:
#            if isinstance(u, pymet.BasicReaction):
#                self.reactions.remove(u)
#            elif isinstance(u, pymet.BasicCompound):
#                self.compounds.remove(u)
#            else:
#                raise TypeError("unidentified metabolic type '{0}'".format(type(u)))
#        if v not in self:
#            if isinstance(v, pymet.BasicReaction):
#                self.reactions.remove(v)
#            elif isinstance(v, pymet.BasicCompound):
#                self.compounds.remove(v)
#            else:
#                raise TypeError("unidentified metabolic type '{0}'".format(type(v)))

    def add_node(self, n, **kw_args):
        """
        """
        if isinstance(n, pymet.BasicReaction):
            self.reactions.add(n)
        elif isinstance(n, pymet.BasicCompound):
            self.compounds.add(n)
        else:
            raise TypeError("unidentified metabolic type '{0}'".format(type(n)))
        super(MetabolicNetwork, self).add_node(n, **kw_args)

    def remove_node(self, n):
        """
        """
        if isinstance(n, pymet.BasicReaction):
            self.reactions.remove(n)
        elif isinstance(n, pymet.BasicCompound):
            self.compounds.remove(n)
        else:
            raise TypeError("unidentified metabolic type '{0}'".format(type(n)))
        super(MetabolicNetwork, self).remove_node(n)

    def remove_compartment(self, compartment):
        rm = set()
        if hasattr(compartment, "__iter__"):
            for rxn in self.reactions:
                if any(cmpd in cmprtmnt for cmpd in rxn.compounds_iter() for
                        cmprtmnt in compartment):
                    rm.add(rxn)
            for rxn in rm:
                self.remove_node(rxn)
            for cmprtmnt in compartment:
                for cmpd in cmprtmnt:
                    self.remove_node(cmpd)
            self.compartments.difference_update(set(compartment))
        else:
            for rxn in self.reactions:
                if any(cmpd in compartment for cmpd in rxn.compounds_iter()):
                    rm.add(rxn)
            for rxn in rm:
                self.remove_node(rxn)
            for cmpd in compartment:
                    self.remove_node(cmpd)
            self.compartments.remove(compartment)

    def to_pruned_currency(self, kegg_information, rpair_categories=["main"],
            scl_threshold=0.4):
        """
        A method to remove undesirable links in metabolic networks based on the
        combination of two different approaches.

        This method applies the "strength of chemical linkage" (SCL) concept [1]_
        to a network seeded with reaction partner information [2]_. Since the
        goal is not to search the shortest biochemically sensible path between
        two compounds but rather the pruning of undesired links involving
        currency metabolites, we use the KEGG RPAIR information, where
        available, to identify relevant interactions between partners in
        reactions. Unavailable KEGG information is supplemented with the SCL
        concept to identify candidates for interactions.

        References
        ----------
        .. [1] Zhou, W., Nakhleh, L., 2011.
           "The strength of chemical linkage as a criterion for pruning
           metabolic graphs".
           Bioinformatics 27, 1957 –1963.

        .. [2] Faust, K., Croes, D., van Helden, J., 2009.
           "Metabolic Pathfinding Using RPAIR Annotation".
           Journal of Molecular Biology 388, 390–414.

        """

        def scl(cmpd_a, cmpd_d):
            scl = 0
            a_total = 0
            d_total = 0
            keys = set(cmpd_a.formula).union(set(cmpd_d.formula))
            if "H" in keys:
                keys.remove("H")
            for atom in keys:
                a_freq = cmpd_a.formula.get(atom, 0)
                d_freq = cmpd_d.formula.get(atom, 0)
                a_total += a_freq
                d_total += d_freq
                scl += min(a_freq, d_freq)
            return float(scl) / float(max(a_total, d_total))

        no_info = set()
        connections = dict()
        # first part, seed the network with kegg rpair information
        for rxn in self.reactions:
            ec_number = rxn.notes.get("ec_number", False)
            if not ec_number:
                no_info.add(rxn)
                continue
            rpairs = False
            for kegg_rxn in kegg_information:
                if kegg_rxn.enzyme == ec_number:
                    rpairs = kegg_rxn.rpair
                    break
            if not rpairs:
                no_info.add(rxn)
                continue
            connections[rxn] = list()
            for cat in rpair_categories:
                pairs = rpairs.get(cat, [])
                for (u, v) in pairs:
                    for cmpd in self.compounds:
                        if cmpd.kegg_id == u and self.has_edge(cmpd, rxn):
                            u = cmpd
                            break
                    for cmpd in self.compounds:
                        if cmpd.kegg_id == v and self.has_edge(rxn, cmpd):
                            v = cmpd
                            break
                    if isinstance(u, pymet.BasicCompound) and isinstance(v,
                            pymet.BasicCompound):
                        connections[rxn].append((u, v))
            if len(connections[rxn]) == 0:
                connections.pop(rxn)
                no_info.add(rxn)
        length = len(no_info)
        last = length + 1
        # find reactions with info just one hop from those without
        while (length - last) < 0 and length > 0:
            for rxn in no_info:
                connections[rxn] = list()
                sources = defaultdict(set)
                for cmpd in self.predecessors_iter(rxn):
                    for src_rxn in self.predecessors_iter(cmpd):
                        if src_rxn in connections:
                            if cmpd in [pair[1] for pair in connections[src_rxn]]:
                                sources[cmpd].add(src_rxn)
                                continue
                            if src_rxn.reversible and cmpd in\
                                    [pair[0] for pair in connections[src_rxn]]:
                                sources[cmpd].add(src_rxn)
                                continue
                for cmpd in sources:
                    for prod in self.successors_iter(rxn):
                        connection = connections.get(rxn, list())
                        if (cmpd, prod) in connection:
                            continue
                        val = scl(cmpd, prod)
                        if val >= scl_threshold:
                            connection.append((cmpd, prod))
                            connections[rxn] = connection
                targets = defaultdict(set)
                for cmpd in self.successors_iter(rxn):
                    for tar_rxn in self.successors_iter(cmpd):
                        if tar_rxn in connections:
                            if cmpd in [pair[0] for pair in connections[tar_rxn]]:
                                targets[cmpd].add(tar_rxn)
                                continue
                            if tar_rxn.reversible and cmpd in\
                                    [pair[1] for pair in connections[tar_rxn]]:
                                targets[cmpd].add(tar_rxn)
                                continue
                for cmpd in targets:
                    for subs in self.predecessors_iter(rxn):
                        connection = connections.get(rxn, list())
                        if (subs, cmpd) in connection:
                            continue
                        val = scl(subs, cmpd)
                        if val >= scl_threshold:
                            connection.append((subs, cmpd))
                            connections[rxn] = connection
                if not connections[rxn]:
                    connections.pop(rxn)
            no_info.difference_update(set(connections))
            last = length
            length = len(no_info)
        if len(connections) != len(self.reactions):
            LOGGER.warn("insufficient reactant pair information, difference of %d",
                    abs(len(connections) - len(self.reactions)))
        pruned = MetabolicNetwork()
        pruned.graph = self.graph.copy()
        pruned.name = "pruned_currency_" + self.name
        pruned.compartments = self.compartments.copy()
        for cmpd in self.compounds:
            pruned.add_node(cmpd)
        for rxn in self.reactions:
            pruned.add_node(rxn)
        for (rxn, rpairs) in connections.iteritems():
            for (u, v) in rpairs:
                pruned.add_edge(u, rxn, **self.edge[u][rxn].copy())
                pruned.add_edge(rxn, v, **self.edge[rxn][v].copy())
        pruned.rpairs = connections
        return pruned

    def to_decompartmentalized(self):
        net = MetabolicNetwork(name="decompartmentalized " + self.name,
                **self.graph)
        for rxn in self.reactions:
            # substrates as reaction attribute
            members = dict()
            for (cmpd, factor) in rxn.substrates.iteritems():
                if isinstance(cmpd, pymet.BasicCompartmentCompound):
                    members[cmpd.compound] = factor
                else:
                    members[cmpd] = factor
            rxn.substrates = members
            # substrates in topological sense
            for cmpd in self.predecessors_iter():
                if isinstance(cmpd, pymet.BasicCompartmentCompound):
                    net.add_edge(cmpd.compound, rxn, **self.edge[cmpd][rxn].copy())
                else:
                    net.add_edge(cmpd, rxn, **self.edge[cmpd][rxn].copy())
            # products as reaction attribute
            members = dict()
            for (cmpd, factor) in rxn.products.iteritems():
                if isinstance(cmpd, pymet.BasicCompartmentCompound):
                    members[cmpd.compound] = factor
                else:
                    members[cmpd] = factor
            rxn.products = members
            # products in topological sense
            for cmpd in self.successors_iter():
                if isinstance(cmpd, pymet.BasicCompartmentCompound):
                    net.add_edge(rxn, cmpd.compound, **self.edge[cmpd][rxn].copy())
                else:
                    net.add_edge(rxn, cmpd, **self.edge[cmpd][rxn].copy())
        self.compartments = set()
        return net

    def introduce_bidirectional(self):

        def add_rev(u, v):
            template.add_edge(u, v)
            template.add_edge(v, u)

        def add_single(u, v):
            template.add_edge(u, v)

        template = MetabolicNetwork(name=self.name)
        for cmpd in self.compounds:
            template.add_node(cmpd)
        for rxn in self.reactions:
            template.add_node(rxn)
        # introduce bidirectional edges
        for rxn in self.reactions:
            if rxn.reversible:
                new_edge = add_rev
            else:
                new_edge = add_single
            for cmpd in self.predecessors_iter(rxn):
                new_edge(cmpd, rxn)
            for cmpd in self.successors_iter(rxn):
                new_edge(rxn, cmpd)
        return template

    def to_compound_centric(self, bidirectional=True):
        """
        """

        def add_bi(u, v, **attr):
            network.add_edge(u, v, **attr)
            network.add_edge(v, u, **attr)

        network = CompoundCentricMultiNetwork()
        network.graph = self.graph.copy()
        network.name="compound_centric_" + self.name
        network.compartments = self.compartments.copy()
        # project to unipartite network with only compound nodes
        for cmpd in self.compounds:
            network.add_node(cmpd)
        # if available, we only add reactant pairs to the network
        if self.rpairs:
            for (rxn, rpairs) in self.rpairs.iteritems():
                if bidirectional and rxn.reversible:
                    # add a bidirectional link
                    add_link = add_bi
                else:
                    # add a unidirectional link
                    add_link = network.add_edge
                for (u, v) in rpairs:
                    add_link(u, v, key=rxn)
        else:
            for rxn in self.reactions:
                if bidirectional and rxn.reversible:
                    # add a bidirectional link
                    add_link = add_bi
                else:
                    # add a unidirectional link
                    add_link = network.add_edge
                for pred in self.predecessors_iter(rxn):
                    for succ in self.successors_iter(rxn):
                        add_link(pred, succ, key=rxn)
        network.remove_edges_from(network.selfloop_edges())
        return network

    def to_reaction_centric(self, bidirectional=True):
        """
        """
        network = ReactionCentricMultiNetwork()
        network.graph = self.graph.copy()
        network.name="reaction_centric_" + self.name
        network.compartments = self.compartments.copy()
        # project to unipartite network with only reaction nodes
        for rxn in self.reactions:
            network.add_node(rxn)
        if self.rpairs:
            for (rxn, rpairs) in self.rpairs.iteritems():
                for (u, v) in rpairs:
        # check whether u is in a reactant pair of a preceeding reaction
                    for src_rxn in self.predecessors_iter(u):
                        if u in [pair[1] for pair in self.rpairs[src_rxn]]:
                            network.add_edge(src_rxn, rxn, key=u)
                        if src_rxn.reversible and u in\
                                [pair[0] for pair in self.rpairs[src_rxn]]:
                            network.add_edge(src_rxn, rxn, key=u)
                        if rxn.reversible and src_rxn.reversible and\
                                network.has_edge(src_rxn, rxn, key=u):
                            network.add_edge(rxn, src_rxn, key=u)
        # check whether u is in a reactant pair of a following reversible reaction
                    for src_rxn in self.successors_iter(u):
                        if src_rxn == rxn:
                            continue
                        if src_rxn.reversible and u in\
                                [pair[0] for pair in self.rpairs[src_rxn]]:
                            network.add_edge(src_rxn, rxn, key=u)
                        if rxn.reversible and src_rxn.reversible and\
                                network.has_edge(src_rxn, rxn, key=u):
                            network.add_edge(rxn, src_rxn, key=u)
        # check whether v is in a reactant pair of a following reaction
                    for tar_rxn in self.successors_iter(v):
                        if v in [pair[0] for pair in self.rpairs[tar_rxn]]:
                            network.add_edge(rxn, tar_rxn, key=v)
                        if tar_rxn.reversible and v in\
                                [pair[1] for pair in self.rpairs[tar_rxn]]:
                            network.add_edge(rxn, tar_rxn, key=v)
                        if rxn.reversible and tar_rxn.reversible and\
                                network.has_edge(rxn, tar_rxn, key=v):
                            network.add_edge(tar_rxn, rxn, key=v)
        # check whether v is in a reactant pair of a preceeding reversible reaction
                    for tar_rxn in self.predecessors_iter(v):
                        if tar_rxn == rxn:
                            continue
                        if tar_rxn.reversible and v in\
                                [pair[1] for pair in self.rpairs[tar_rxn]]:
                            network.add_edge(rxn, tar_rxn, key=v)
                        if rxn.reversible and tar_rxn.reversible and\
                                network.has_edge(rxn, tar_rxn, key=v):
                            network.add_edge(tar_rxn, rxn, key=v)
        else:
            for cmpd in self.compounds:
                predecessors = self.predecessors(cmpd)
                successors = self.successors(cmpd)
                if bidirectional:
                    rev_pred = [rxn for rxn in predecessors if rxn.reversible]
                    rev_succ = [rxn for rxn in successors if rxn.reversible]
                for pred in predecessors:
                    for succ in successors:
                        network.add_edge(pred, succ, key=cmpd)
                # add links due to reversibility
                    for rxn in rev_pred:
                        network.add_edge(pred, rxn, key=cmpd)
                for rxn in rev_succ:
                    for succ in successors:
                        network.add_edge(rxn, succ, key=cmpd)
                    for pred in rev_pred:
                        network.add_edge(rxn, pred, key=cmpd)
        # we added a lot of self-links in the process, I felt removing them
        # later was more efficient than working with set differences all the
        # time
        network.remove_edges_from(network.selfloop_edges())
        return network

    def draw(self, filename, output_format="pdf", layout_program="fdp",
                layout_args="", distinct=False):
        import pygraphviz as pgv
        OPTIONS = misc.OptionsManager.get_instance()
        net = pgv.AGraph(directed=True, name=filename, strict=True)
        node_attr= dict()
        link_attr= dict()
        # add compound nodes
        indeces = dict(itertools.izip(self.compounds, itertools.count()))
        for (cmpd, i) in indeces.iteritems():
            net.add_node(i, label=str(cmpd), shape="ellipse", **node_attr)
        # add reactions
        indeces.update(itertools.izip(self.reactions,
                itertools.count(len(self.compounds))))
        i = len(self.compounds) + len(self.reactions)
        for rxn in self.reactions:
            net.add_node(indeces[rxn], label=str(rxn), shape="box", **node_attr)
            # add forward reaction links
            for cmpd in self.predecessors(rxn):
                net.add_edge(indeces[cmpd], indeces[rxn], **link_attr)
            for cmpd in self.successors(rxn):
                net.add_edge(indeces[rxn], indeces[cmpd], **link_attr)
            if rxn.reversible:
                if distinct:
                    rev = pymet.BasicReaction(str(rxn) + OPTIONS.reversible_suffix)
                    indeces[rev] = i
                    net.add_node(i, label=str(rev), shape="box", **node_attr)
                    # add backward reaction links
                    for cmpd in self.predecessors(rxn):
                        net.add_edge(indeces[rev], indeces[cmpd], **link_attr)
                    for cmpd in self.successors(rxn):
                        net.add_edge(indeces[cmpd], indeces[rev], **link_attr)
                    i += 1
                else:
                    # add backward reaction links
                    for cmpd in self.predecessors(rxn):
                        net.add_edge(indeces[rxn], indeces[cmpd],
                                style="dotted", **link_attr)
                    for cmpd in self.successors(rxn):
                        net.add_edge(indeces[cmpd], indeces[rxn],
                                style="dotted", **link_attr)
        filename = "%s.%s" % (filename, output_format)
        net.draw(filename, prog=layout_program, args=layout_args)

    def to_system(self):
        system = pymet.MetabolicSystem(name=self.name)
        for rxn in self.reactions:
            subs = dict((pymet.SBMLCompound(str(cmpd)),
                    self[cmpd][rxn]["coefficient"]) for cmpd in self.pred[rxn])
            prods = dict((pymet.SBMLCompound(str(cmpd)),
                    self[rxn][cmpd]["coefficient"]) for cmpd in self.succ[rxn])
            if rxn.reversible:
                system.add(pymet.SBMLReaction(str(rxn), subs, prods,
                        rxn.reversible, lower_bound=-OPTIONS.upper_bound,
                        upper_bound=OPTIONS.upper_bound))
            else:
                system.add(pymet.SBMLReaction(str(rxn), subs, prods,
                        rxn.reversible, lower_bound=OPTIONS.lower_bound,
                        upper_bound=OPTIONS.upper_bound))
        return system


class CompoundCentricNetwork(nx.DiGraph):
    """
    """

    def __init__(self, data=None, name="", compartments=None, **kw_args):
        """
        """
        super(CompoundCentricNetwork, self).__init__(data=data, name=name, **kw_args)
        self.compartments = misc.convert(compartments, set, set())

    def draw(self, filename, output_format="pdf", layout_program="fdp", layout_args=""):
        import pygraphviz as pgv
        net = pgv.AGraph(directed=True, name=filename, strict=True)
        node_attr= dict()
        link_attr= dict()
        indeces = dict()
        # add compound nodes
        for (i, cmpd) in enumerate(self.nodes_iter()):
            indeces[cmpd] = i
            net.add_node(i, label=str(cmpd), shape="ellipse", **node_attr)
        # add links
        for (u, v) in self.edges_iter():
            net.add_edge(indeces[u], indeces[v], **link_attr)
        filename = "%s.%s" % (filename, output_format)
        net.draw(filename, prog=layout_program, args=layout_args)


class CompoundCentricMultiNetwork(nx.MultiDiGraph):
    """
    """

    def __init__(self, data=None, name="", compartments=None, **kw_args):
        """
        """
        super(CompoundCentricMultiNetwork, self).__init__(data=data, name=name, **kw_args)
        self.compartments = misc.convert(compartments, set, set())

    def to_directed(self):
        """
        Return a copy with no multiple edges and no attributes.
        """
        network = CompoundCentricNetwork()
        network.graph = self.graph.copy()
        network.name="directed_" + self.name
        network.compartments = self.compartments.copy()
        network.add_nodes_from(self.nodes_iter())
        network.add_edges_from(self.edges_iter())
        return network

    def draw(self, filename, output_format="pdf", layout_program="fdp", layout_args=""):
        import pygraphviz as pgv
        net = pgv.AGraph(directed=True, name=filename, strict=False)
        node_attr= dict()
        link_attr= dict()
        indeces = dict()
        # add compound nodes
        for (i, cmpd) in enumerate(self.nodes_iter()):
            indeces[cmpd] = i
            net.add_node(i, label=str(cmpd), shape="ellipse", **node_attr)
        # add links
        for (u, v) in self.edges_iter():
            net.add_edge(indeces[u], indeces[v], **link_attr)
        filename = "%s.%s" % (filename, output_format)
        net.draw(filename, prog=layout_program, args=layout_args)


class ReactionCentricNetwork(nx.DiGraph):
    """
    """

    def __init__(self, data=None, name="", compartments=None, **kw_args):
        """
        """
        super(ReactionCentricNetwork, self).__init__(data=data, name=name, **kw_args)
        self.compartments = misc.convert(compartments, set, set())

    def draw(self, filename, output_format="pdf", layout_program="fdp", layout_args=""):
        """
        """
        import pygraphviz as pgv
        net = pgv.AGraph(directed=True, name=filename, strict=False)
        node_attr= dict()
        link_attr= dict()
        indeces = dict()
        # add reaction nodes
        for (i, rxn) in enumerate(self.nodes_iter()):
            indeces[rxn] = i
            net.add_node(i, label=str(rxn), shape="box", **node_attr)
        # add links
        for (u, v) in self.edges_iter():
            net.add_edge(indeces[u], indeces[v], **link_attr)
        filename = "%s.%s" % (filename, output_format)
        net.draw(filename, prog=layout_program, args=layout_args)


class ReactionCentricMultiNetwork(nx.MultiDiGraph):
    """
    """

    def __init__(self, data=None, name="", compartments=None, **kw_args):
        """
        """
        super(ReactionCentricMultiNetwork, self).__init__(data=data, name=name, **kw_args)
        self.compartments = misc.convert(compartments, set, set())

    def to_directed(self):
        """
        Return a copy with no multiple edges and no attributes.
        """
        network = ReactionCentricNetwork()
        network.graph = self.graph.copy()
        network.name="directed_" + self.name
        network.compartments = self.compartments.copy()
        network.add_nodes_from(self.nodes_iter())
        network.add_edges_from(self.edges_iter())
        return network

    def draw(self, filename, output_format="pdf", layout_program="fdp", layout_args=""):
        import pygraphviz as pgv
        net = pgv.AGraph(directed=True, name=filename, strict=False)
        node_attr= dict()
        link_attr= dict()
        indeces = dict()
        # add reaction nodes
        for (i, rxn) in enumerate(self.nodes_iter()):
            indeces[rxn] = i
            net.add_node(i, label=str(rxn), shape="box", **node_attr)
        # add links
        for (u, v) in self.edges_iter():
            net.add_edge(indeces[u], indeces[v], **link_attr)
        filename = "%s.%s" % (filename, output_format)
        net.draw(filename, prog=layout_program, args=layout_args)


def read_edgelist(filename, sep="\t", comment="#", mode="rb",
        encoding="utf-8", **kw_args):
    """
    """
    def build_node(name):
        if name.startswith(OPTIONS.compound_prefix):
            compound = pymet.BasicCompound(name[len(OPTIONS.compound_prefix):])
            return compound
        elif name.startswith(OPTIONS.reaction_prefix):
            if name.endswith(OPTIONS.reversible_suffix):
                reaction = pymet.BasicReaction(name[len(OPTIONS.reaction_prefix):
                        -len(OPTIONS.reversible_suffix)])
                reaction.reversible = True
            else:
                reaction = pymet.BasicReaction(name[len(OPTIONS.reaction_prefix):])
            return reaction
        else:
            raise PyOrganismError("unidentified metabolic object '{0}'".format(name))

    net = MetabolicNetwork(name=filename)
    # read the file contents
    kw_args["mode"] = mode
    kw_args["encoding"] = encoding
    with  open_file(filename, **kw_args) as (file_h, ext):
        # parse the edgelist into a simple metabolic network
        for row in read_tabular(file_h, sep=sep, comment=comment):
            u = build_node(row[0])
            if isinstance(u, pymet.BasicReaction) and\
                    row[0].endswith(OPTIONS.reversible_suffix):
                continue
            v = build_node(row[1])
            if isinstance(v, pymet.BasicReaction) and\
                    row[1].endswith(OPTIONS.reversible_suffix):
                continue
            net.add_edge(u, v)
    return net

def write_edgelist(network, filename, distinct=True, delimiter="\t",
        comment="#", mode="wb", encoding="utf-8", **kw_args):
    """
    """
    lines = list()
    # assemble lines
    for rxn in network.reactions:
        rxn_name = OPTIONS.reaction_prefix + rxn.name
        if rxn.reversible:
            if distinct:
                rev_name = "{0}{1}{2}".format(OPTIONS.reaction_prefix, rxn.name,
                        OPTIONS.reversible_suffix)
            else:
                rev_name = rxn_name
            for cmpd in network.successors_iter(rxn):
                lines.append("{0}{1}{2}\n".format(rxn_name, delimiter,
                        OPTIONS.compound_prefix + cmpd.name))
                lines.append("{0}{1}{2}\n".format(OPTIONS.compound_prefix + cmpd.name,
                        delimiter, rev_name))
            for cmpd in network.predecessors_iter(rxn):
                lines.append("{0}{1}{2}\n".format(OPTIONS.compound_prefix + cmpd.name,
                        delimiter, rxn_name))
                lines.append("{0}{1}{2}\n".format(rev_name, delimiter,
                        OPTIONS.compound_prefix + cmpd.name))
        else:
            for cmpd in network.successors_iter(rxn):
                lines.append("{0}{1}{2}\n".format(rxn_name, delimiter,
                        OPTIONS.compound_prefix + cmpd.name))
            for cmpd in network.predecessors_iter(rxn):
                lines.append("{0}{1}{2}\n".format(OPTIONS.compound_prefix + cmpd.name,
                        delimiter, rxn_name))
    # write to file
    kw_args["mode"] = mode
    kw_args["encoding"] = encoding
    with  open_file(filename, **kw_args) as (file_h, ext):
        file_h.writelines(lines)

