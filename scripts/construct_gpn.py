# -*- coding: utf-8 -*-


import sys
import os
import logging

import networkx as nx

import pyorganism
import pyorganism.regulation as pyreg


LOGGER = logging.getLogger()
LOGGER.addHandler(logging.StreamHandler())
LOGGER.setLevel(logging.INFO)


def construct_gpn(path, window_sizes):
    LOGGER.info("{0:*^78s}".format("Construct GPN"))
    version = os.path.basename(path)
    if not version:
        version = os.path.basename(os.path.dirname(path))
    LOGGER.info("{0:*^78s}".format(version))
    genes = pyorganism.read_pickle(os.path.join(path, "genes.pkl"))
    gpn_gen = pyreg.GPNGenerator()
    gpn_gen.parse_genes(genes)
    for window in window_sizes:
        LOGGER.info("Window Size: %d",window)
        gpn = gpn_gen.generate_gpn(window)
        LOGGER.info("%d Nodes", len(gpn))
        LOGGER.info("%d Links", gpn.size())
        components = list(nx.connected_components(gpn))
        LOGGER.info("%d Components", len(components))
        LOGGER.info("Largest Component: %d", len(components[0]))
        if len(components) > 1:
            LOGGER.info("Second Largest Component: %d", len(components[1]))
        pyorganism.write_pickle(gpn, os.path.join(path, "gpn_%d.pkl" % window))

if __name__ == "__main__":
    if len(sys.argv) != 2:
        LOGGER.critical("%s <RegulonDB objects path>", sys.argv[0])
        sys.exit(2)
    else:
        construct_gpn(sys.argv[1], [4000, 5000, 6000, 7000])

