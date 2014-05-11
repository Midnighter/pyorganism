# -*- coding: utf-8 -*-


import sys
import os
import logging

import numpy
import networkx as nx

import pyorganism
import pyorganism.regulation as pyreg


LOGGER = logging.getLogger()
LOGGER.addHandler(logging.StreamHandler())
LOGGER.setLevel(logging.INFO)


class NodeConverter(object):
    def __init__(self, t_units, default="", **kw_args):
        """
        Sets up an instance that makes one attribute of a collection of objects
        searcheable.

        These attributes can have any value that is comparable. In none string
        cases, a different default value should be passed.
        """
        super(NodeConverter, self).__init__(**kw_args)
        self.t_units = list(t_units)
        self.targets = numpy.array([tu.promoter.unique_id if tu.promoter else default\
                for tu in self.t_units])
        self.indeces = numpy.arange(len(self.targets), dtype=int)

    def __call__(self, value):
        return (gene for i in self.indeces[self.targets == value] for gene in self.t_units[i].genes)


def construct_trn(path):
    version = os.path.basename(path)
    if not version:
        version = os.path.basename(os.path.dirname(path))
    LOGGER.info("{0:*^78s}".format(version))
    # load objects so that they are in memory
    conformations = pyorganism.read_pickle(os.path.join(path, "conformations.pkl"))
    t_units = pyorganism.read_pickle(os.path.join(path, "transcription_units.pkl"))
    interactions = pyorganism.read_pickle(os.path.join(path, "interactions.pkl"))
    assert all(triple[0] in pyreg.Conformation for triple in interactions),\
            "unknown conformation in regulatory interactions"
    assert all(triple[1] in pyreg.Promoter for triple in interactions),\
            "unknown promoter in regulatory interactions"
    # conformation-promoter directed regulatory network
    cpn = nx.MultiDiGraph()
    for (u, v, inter) in interactions:
        first = pyreg.Conformation[u]
        second = pyreg.Promoter[v]
        cpn.add_edge(first, second, key=inter)
    # TRN from CPN
    trn = pyreg.TRN()
    node_converter = NodeConverter(t_units)
    for (u, v, k, d) in cpn.edges_iter(keys=True, data=True):
        first = u.t_factor
        for nbr in node_converter(v.unique_id):
            trn.add_edge(first, nbr, key=k, **d.copy())
    LOGGER.info("%d Nodes", len(trn))
    LOGGER.info("%d Links", trn.size())
    components = nx.connected_components(trn.to_undirected())
    LOGGER.info("%d Components", len(components))
    LOGGER.info("Largest Component: %d", len(components[0]))
    if len(components) > 1:
        LOGGER.info("Second Largest Component: %d", len(components[1]))
    pyorganism.write_pickle(trn, os.path.join(path, "trn.pkl"))


if __name__ == "__main__":
    if len(sys.argv) != 2:
        LOGGER.critical("%s <RegulonDB objects path>", sys.argv[0])
        sys.exit(2)
    else:
        construct_trn(sys.argv[1])

