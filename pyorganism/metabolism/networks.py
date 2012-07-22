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


__all__ = ["CompoundCentricNetwork", "CompoundCentricMultiNetwork",
        "ReactionCentricNetwork", "ReactionCentricMultiNetwork",
        "MetabolicNetwork"]

import logging
import itertools
import networkx as nx

from . import elements as pymet
from .. import miscellaneous as misc


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
        self.compartments = misc.convert(compartments, set)
        self.compounds = set()
        self.reactions = set()

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

    def to_compound_centric(self):
        """
        """

        def add_bi(u, v, **attr):
            network.add_edge(u, v, **attr)
            network.add_edge(v, u, **attr)

        network = CompoundCentricMultiNetwork(name="compound_centric_" + self.name)
        # project to unipartite network with only compound nodes
        for cmpd in self.compounds:
            network.add_node(cmpd)
        for rxn in self.reactions:
            if rxn.reversible:
                # add a bidirectional link
                add_link = add_bi
            else:
                # add a unidirectional link
                add_link = network.add_edge
            for pred in self.predecessors_iter(rxn):
                for succ in self.successors_iter(rxn):
                    add_link(pred, succ, reaction=rxn)
        network.remove_edges_from(network.selfloop_edges())
        return network

    def to_reaction_centric(self):
        """
        """
        network = ReactionCentricMultiNetwork(name="reaction_centric_" + self.name)
        # project to unipartite network with only reaction nodes
        for rxn in self.reactions:
            network.add_node(rxn)
        for cmpd in self.compounds:
            predecessors = self.predecessors(cmpd)
            rev_pred = [rxn for rxn in predecessors if rxn.reversible]
            successors = self.successors(cmpd)
            rev_succ = [rxn for rxn in successors if rxn.reversible]
            for pred in predecessors:
                for succ in successors:
                    network.add_edge(pred, succ, compound=cmpd)
            # add links due to reversibility
                for rxn in rev_pred:
                    network.add_edge(pred, rxn, compound=cmpd)
            for rxn in rev_succ:
                for succ in successors:
                    network.add_edge(rxn, succ, compound=cmpd)
                for pred in rev_pred:
                    network.add_edge(rxn, pred, compound=cmpd)
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
            net.add_node(i, label=cmpd.name, shape="circle", **node_attr)
        # add reactions
        indeces.update(itertools.izip(self.reactions,
                itertools.count(len(self.compounds))))
        i = len(self.compounds) + len(self.reactions)
        for rxn in self.reactions:
            net.add_node(indeces[rxn], label=rxn.name, shape="box", **node_attr)
            # add forward reaction links
            for cmpd in self.predecessors(rxn):
                net.add_edge(indeces[cmpd], indeces[rxn], **link_attr)
            for cmpd in self.successors(rxn):
                net.add_edge(indeces[rxn], indeces[cmpd], **link_attr)
            if rxn.reversible:
                if distinct:
                    rev = pymet.BasicReaction(rxn.name + OPTIONS.reversible_suffix)
                    indeces[rev] = i
                    net.add_node(i, label=rev.name, shape="box", **node_attr)
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

    def __init__(self, data=None, name="", **kw_args):
        """
        """
        super(CompoundCentricNetwork, self).__init__(data=data, name=name, **kw_args)

    def draw(self, filename, output_format="pdf", layout_program="fdp", layout_args=""):
        import pygraphviz as pgv
        net = pgv.AGraph(directed=True, name=filename, strict=True)
        node_attr= dict()
        link_attr= dict()
        indeces = dict()
        # add compound nodes
        for (i, cmpd) in enumerate(self.nodes_iter()):
            indeces[cmpd] = i
            net.add_node(i, label=cmpd.name, shape="circle", **node_attr)
        # add links
        for (u, v) in self.edges_iter():
            net.add_edge(indeces[u], indeces[v], **link_attr)
        filename = "%s.%s" % (filename, output_format)
        net.draw(filename, prog=layout_program, args=layout_args)


class CompoundCentricMultiNetwork(nx.MultiDiGraph):
    """
    """

    def __init__(self, data=None, name="", **kw_args):
        """
        """
        super(CompoundCentricMultiNetwork, self).__init__(data=data, name=name, **kw_args)

    def to_directed(self):
        """
        Return a copy with no multiple edges and no attributes.
        """
        copy = CompoundCentricNetwork(name="directed " + self.name)
        copy.add_nodes_from(self.nodes_iter())
        copy.add_edges_from(self.edges_iter())
        return copy

    def draw(self, filename, output_format="pdf", layout_program="fdp", layout_args=""):
        import pygraphviz as pgv
        net = pgv.AGraph(directed=True, name=filename, strict=False)
        node_attr= dict()
        link_attr= dict()
        indeces = dict()
        # add compound nodes
        for (i, cmpd) in enumerate(self.nodes_iter()):
            indeces[cmpd] = i
            net.add_node(i, label=cmpd.name, shape="circle", **node_attr)
        # add links
        for (u, v) in self.edges_iter():
            net.add_edge(indeces[u], indeces[v], **link_attr)
        filename = "%s.%s" % (filename, output_format)
        net.draw(filename, prog=layout_program, args=layout_args)


class ReactionCentricNetwork(nx.DiGraph):
    """
    """

    def __init__(self, data=None, name="", **kw_args):
        """
        """
        super(ReactionCentricNetwork, self).__init__(data=data, name=name, **kw_args)

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
            net.add_node(i, label=rxn.name, shape="box", **node_attr)
        # add links
        for (u, v) in self.edges_iter():
            net.add_edge(indeces[u], indeces[v], **link_attr)
        filename = "%s.%s" % (filename, output_format)
        net.draw(filename, prog=layout_program, args=layout_args)


class ReactionCentricMultiNetwork(nx.MultiDiGraph):
    """
    """

    def __init__(self, data=None, name="", **kw_args):
        """
        """
        super(ReactionCentricMultiNetwork, self).__init__(data=data, name=name, **kw_args)

    def to_directed(self):
        """
        Return a copy with no multiple edges and no attributes.
        """
        copy = CompoundCentricNetwork(name="directed_" + self.name)
        copy.add_nodes_from(self.nodes_iter())
        copy.add_edges_from(self.edges_iter())
        return copy

    def draw(self, filename, output_format="pdf", layout_program="fdp", layout_args=""):
        import pygraphviz as pgv
        net = pgv.AGraph(directed=True, name=filename, strict=False)
        node_attr= dict()
        link_attr= dict()
        indeces = dict()
        # add reaction nodes
        for (i, rxn) in enumerate(self.nodes_iter()):
            indeces[rxn] = i
            net.add_node(i, label=rxn.name, shape="box", **node_attr)
        # add links
        for (u, v) in self.edges_iter():
            net.add_edge(indeces[u], indeces[v], **link_attr)
        filename = "%s.%s" % (filename, output_format)
        net.draw(filename, prog=layout_program, args=layout_args)

