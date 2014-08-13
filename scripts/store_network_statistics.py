#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
=====================================
Unified Storage of Network Statistics
=====================================

:Authors:
    Moritz Emanuel Beber
:Date:
    2014-01-29
:Copyright:
    Copyright(c) 2014 Jacobs University of Bremen. All rights reserved.
:File:
    store_network_statistics.py
"""


# stdlib
import sys
import os
import logging
import itertools
# external
import pandas
import networkx as nx
# project
import pyorganism
import pyorganism.regulation as pyreg # needed for object definitions

from pyorganism.io.regulondb import RELEASE
from meb.utils.network.subgraphs import triadic_census


LOGGER = logging.getLogger()
LOGGER.addHandler(logging.StreamHandler())
LOGGER.setLevel(logging.INFO)

def trn_stats(genes, trn, version):
    LOGGER.info("Computing TRN statistics")
    grn = trn.to_grn()
    nodes = sorted(grn.nodes())
    regulating = set(node for (node, deg) in grn.out_degree_iter() if deg > 0)
    regulated = set(node for (node, deg) in grn.in_degree_iter() if deg > 0)
    components = list(nx.weakly_connected_components(grn))
    data = dict()
    for (a, b) in itertools.product(("in", "out"), repeat=2):
        data["{a}_{b}_ass".format(a=a, b=b)] = nx.degree_assortativity_coefficient(grn, x=a, y=b)
    census = triadic_census(grn)
    forward = census["030T"]
    feedback = census["030C"]
    cycles = list(nx.simple_cycles(grn))
    in_deg = [grn.in_degree(node) for node in nodes]
    out_deg = [grn.out_degree(node) for node in nodes]
    bc = nx.betweenness_centrality(grn)
    bc = [bc[node] for node in nodes]
    data["version"] = version
    data["release"] = RELEASE[version]
    data["num_genes"] = len(genes)
    data["num_regulating"] = len(regulating)
    data["num_regulated"] = len(regulated)
    data["num_links"] = grn.size()
    data["num_components"] = len(components)
    data["largest_component"] = len(components[0])
    data["feed_forward"] = forward
    data["feedback"] = feedback
    data["cycles"] = len(cycles)
    data["hub_out_deg"] = max(out_deg)
    stats = pandas.DataFrame(data, index=[1])
    stats["release"] = pandas.to_datetime(stats["release"])
    dists = pandas.DataFrame({
        "version": version,
        "release": RELEASE[version],
        "node": [node.unique_id for node in nodes],
        "in_degree": in_deg,
        "out_degree": out_deg,
        "betweenness": bc
    })
    dists["release"] = stats["release"].copy()
    return (stats, dists)

def gpn_stats(genes, gpn, version):
    LOGGER.info("Computing GPN statistics")
    nodes = sorted(gpn.nodes())
    components = list(nx.connected_components(gpn))
    ass = nx.degree_assortativity_coefficient(gpn)
    deg = [gpn.degree(node) for node in nodes]
    stats = pandas.DataFrame(data={
            "version": version,
            "release": RELEASE[version],
            "num_genes": len(genes),
            "num_nodes": len(gpn),
            "num_links": gpn.size(),
            "num_components": len(components),
            "largest_component": len(components[0]),
            "assortativity": ass,
            "hub_deg": max(deg)
            }, index=[1])
    stats["release"] = pandas.to_datetime(stats["release"])
    dists = pandas.DataFrame(data={
            "version": version,
            "release": RELEASE[version],
            "node": [node.unique_id for node in nodes],
            "degree": deg,
            })
    dists["release"] = stats["release"].copy()
    return (stats, dists)

def store_results(filename, df):
    if os.path.exists(filename):
        results = pandas.read_csv(filename, sep=";", header=0, index_col=False)
        results = pandas.concat([results, df], ignore_index=True)
    else:
        results = df
    results.to_csv(filename, sep=";", header=True, index=False)

def main(in_path, out_path, version=""):
    version = os.path.basename(in_path)
    if not version:
        version = os.path.basename(os.path.dirname(in_path))
    LOGGER.info("{0:*^78s}".format(version))
    LOGGER.info("Loading genes")
    genes = pyorganism.read_pickle(os.path.join(in_path, "genes.pkl"))
    LOGGER.info("Loading TRN")
    trn = pyorganism.read_pickle(os.path.join(in_path, "trn.pkl"))
    LOGGER.info("Loading GPN")
    gpn = pyorganism.read_pickle(os.path.join(in_path, "gpn_5000.pkl"))
    version = os.path.basename(in_path)
    if not version:
        version = os.path.basename(os.path.dirname(in_path))
    (stats, dists) = gpn_stats(genes, gpn, version)
    store_results(os.path.join(out_path, "gpn_5000_statistics.csv"), stats)
    store_results(os.path.join(out_path, "gpn_5000_distributions.csv"), dists)
    (stats, dists) = trn_stats(genes, trn, version)
    store_results(os.path.join(out_path, "trn_statistics.csv"), stats)
    store_results(os.path.join(out_path, "trn_distributions.csv"), dists)


if __name__ == "__main__":
    if len(sys.argv) < 3 or len(sys.argv) > 4:
        LOGGER.critical("%s <RegulonDB objects: path> <Results dir: path>"\
                " [<RegulonDB version: str>]", sys.argv[0])
        sys.exit(2)
    else:
        main(*sys.argv[1:])

