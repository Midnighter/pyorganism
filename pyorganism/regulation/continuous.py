# -*- coding: utf-8 -*-


from __future__ import (absolute_import, unicode_literals, division)


"""
===========================
Regulatory Control Measures
===========================

:Authors:
    Moritz Emanuel Beber
:Date:
    2014-11-01
:Copyright:
    Copyright(c) 2014 Jacobs University of Bremen. All rights reserved.
:File:
    continuous.py
"""


import logging
import warnings
from itertools import chain

import numpy as np
from builtins import dict

from . import shuffling
from . import continuous_wrapper as con
from .elements import (Gene, Regulator, TranscriptionUnit, Operon)
from .. import miscellaneous as misc


__all__ = [
    "ContinuousControl"
]


LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(misc.NullHandler())


class ContinuousControl(object):
    """
    """
    _sampling = {
        "timeline": shuffling.timeline_sample,
        "fork": shuffling.block_sample,
        "fork-strand": shuffling.block_sample,
    }
    _measures = {
        "absolute": con.abs_control_timeline,
        "difference": con.difference_control_timeline,
        "absolute-difference": con.abs_difference_control_timeline,
        "functional": con.functional_control_timeline,
        "functional-comparison": con.functional_comparison_timeline,
        "delayed-functional": con.delayed_functional_control_timeline
    }

    def __init__(self):
        """
        """
        self.num_nodes = None
        self.num_links = None
        self.nodes = None
        self.node2id = None
        self.sources = None
        self.targets = None
        self.functions = None
        self.expression = None
        self.sample = None
        self.sample_args = None
        self.subnet = None

    def setup(self, net, df, feature2node, sampling, node2feature=None, nodes=None,
            oric=(3923767, 3923998),
            size=4639675):
        """
        Convert a pandas DataFrame to an expression matrix.
        """
        if len(net) < 2 or net.size() < 1:
            raise ValueError("aborting due to small network size")
        try:
            self.sample = self._sampling[sampling]
        except KeyError:
            raise ValueError("'{0}' is an unknown sampling"\
                    " method.".format(sampling))
        if nodes is None:
            nodes = sorted(net.nodes())
        if node2feature is None:
            node2feature = {n: feat for (feat, n) in feature2node.items()}
        data = dict()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            for node in nodes:
                if node not in net:
                    continue
                elif isinstance(node, Gene):
                    try:
                        series = df.loc[node2feature[node]].copy()
                    except KeyError:
                        continue
                elif isinstance(node, Regulator):
                    features = {node2feature[gene] for gene in node.coded_from\
                            if gene in node2feature}
                    # unmapped features are None
                    if None in features:
                        features.remove(None)
                    series = np.nanmean(df.loc[features], axis=0)
                elif isinstance(node, (TranscriptionUnit, Operon)):
                    features = {node2feature[gene] for gene in node.genes\
                            if gene in node2feature}
                    # unmapped features are None
                    if None in features:
                        features.remove(None)
                    series = np.nanmean(df.loc[features], axis=0)
                else:
                    raise ValueError("unknown regulatory network node type %r" % (node,))
                if not np.isnan(series).all():
                    data[node] = series
        active = [n for n in nodes if n in data]
        self.subnet = net.subgraph(active)
        nom = len(active)
        denom = len(nodes)
        LOGGER.info("%d/%d active nodes (%.2f%%)", nom, denom, nom / denom * 100.0)
        expression = np.zeros((df.shape[1], len(active)))
        if sampling == "fork":
            ter = int((np.mean(oric) + size / 2) % size)
            rest = list()
            left = list()
            right = list()
            for (i, node) in enumerate(active):
                expression[:, i] = data[node]
                fork_info(i, node, oric, ter, left, right, rest)
            blocks = [0, len(left), len(right), len(rest)]
            self.sample_args = {"blocks": np.cumsum(blocks, dtype=np.int32)}
            # organize expression into complete blocks
            expression = np.ascontiguousarray(expression[:, left + right + rest])
            # reorder nodes to correspond to expression
            self.nodes = [active[i] for i in chain(left, right, rest)]
        elif sampling == "fork-strand":
            ter = int((np.mean(oric) + size / 2) % size)
            rest = list()
            left_forward = list()
            left_reverse = list()
            right_forward = list()
            right_reverse = list()
            for (i, node) in enumerate(active):
                expression[:, i] = data[node]
                fork_strand_info(i, node, oric, ter, left_forward, left_reverse,
                    right_forward, right_reverse, rest)
            blocks = [0, len(left_forward), len(left_reverse),
                    len(right_forward), len(right_reverse), len(rest)]
            self.sample_args = {"blocks": np.cumsum(blocks, dtype=np.int32)}
            # organize expression into complete blocks
            expression = np.ascontiguousarray(expression[:, left_forward + left_reverse +
                right_forward + right_reverse + rest])
            # reorder nodes to correspond to expression
            self.nodes = [active[i] for i in chain(left_forward,
                left_reverse, right_forward, right_reverse, rest)]
        elif sampling == "timeline":
            for (i, node) in enumerate(active):
                expression[:, i] = data[node]
            self.sample_args = dict()
            self.nodes = active
        expression[~np.isfinite(expression)] = 0.0
        self.expression = expression
        self.node2id = {n: i for (i, n) in enumerate(self.nodes)}

    def from_trn(self, function="regulatory"):
        """
        Prepare data structures for continuous control analysis using a
        transcriptional regulatory network.

        Parameters
        ----------
        function: hashable (optional)
            In case of a DiGraph, the keyword for edge data that describes the
            regulatory function. For a MultiDiGraph the key is used.
        """
        self.num_nodes = len(self.nodes)
        sources = list()
        targets = list()
        functions = list()
        if self.subnet.is_multigraph():
            edge_iter = ((u, v, k) for (u, v, k) in self.subnet.edges_iter(keys=True))
        else:
            edge_iter = ((u, v, data[function]) for (u, v, data) in\
                    self.subnet.edges_iter(data=True))
        for (u, v, func) in edge_iter:
            # filter unknown links
            if func == 0:
                continue
            # insert dual regulation as both activating and inhibiting
            elif func == 2:
                sources.append(self.node2id[u])
                targets.append(self.node2id[v])
                functions.append(1)
                sources.append(self.node2id[u])
                targets.append(self.node2id[v])
                functions.append(-1)
                continue
            sources.append(self.node2id[u])
            targets.append(self.node2id[v])
            functions.append(func)
        self.num_links = len(sources)
        self.sources = np.asarray(sources, dtype=np.int32)
        self.targets = np.asarray(targets, dtype=np.int32)
        self.functions = np.asarray(functions, dtype=np.int32)

    def from_gpn(self):
        """
        Prepare data structures for continuous control analysis using a
        gene proximity network.
        """
        self.num_nodes = len(self.nodes)
        self.num_links = self.subnet.size()
        self.sources = np.zeros(self.num_links, dtype=np.int32)
        self.targets = np.zeros(self.num_links, dtype=np.int32)
        for (i, (u, v)) in enumerate(self.subnet.edges_iter()):
            self.sources[i] = self.node2id[u]
            self.targets[i] = self.node2id[v]

    def series_ctc(self, measure, random_num=1E04, delay=0):
        """
        Compute control type confidence for a series of expression levels.

        The expected format of expression values is a matrix, where the first
        dimension corresponds to the points in the series and the second to the
        nodes in the network. The levels for one node should be normalized
        between 0 and 1 over the whole series. The order of nodes must be the same.

        Parameters
        ----------
        expression: numpy.ndarray
            Two dimensional expression levels in the same order as nodes.
        """
        random_num = int(random_num)
        try:
            func = self._measures[measure]
        except KeyError:
            raise ValueError("'{0}' is an unknown control measure."\
                    " Please try one of:\n{1}".format(measure,
                    "\n* ".join(self._measures.keys())))
        if "comparison" in measure:
            expression = np.diff(self.expression, axis=0) # expression was transposed for C code
        else:
            expression = self.expression
        num_points = self.expression.shape[0]
        num_features = self.expression.shape[1]
        assert num_features == self.num_nodes
        control = np.zeros(num_points, dtype=np.double)
        random = np.zeros((random_num, num_points), dtype=np.double)
        kw_args = {
            "sources": self.sources,
            "targets": self.targets,
            "functions": self.functions,
            "num_links": self.num_links,
            "num_points": num_points,
            "dim": num_features,
            "delay": delay
        }
        # compute original control
        kw_args["expression"] = expression
        kw_args["control"] = control
        func(**kw_args)
        # compute randomized control
        for (i, rnd_series) in enumerate(self.sample(expression, random_num,
            **self.sample_args)):
            kw_args["expression"] = rnd_series
            kw_args["control"] = random[i, :]
            func(**kw_args)
        # make points the first dimension
        random = np.ascontiguousarray(random.T)
        z_scores = control - np.nanmean(random, axis=1)
        z_scores /= np.nanstd(random, ddof=1, axis=1)
        return (z_scores, control, random)


def fork_info(i, node, oric, ter, left, right, rest):
    """
    Regulators are ignored for now.
    """
    if isinstance(node, Gene):
        start = node.position_start
    elif isinstance(node, (TranscriptionUnit, Operon)):
        start = min(gene.position_start for gene in node.genes)
    else:
        rest.append(i)
        return
    if start > ter and start < oric[0]:
        right.append(i)
    elif start < ter or start > oric[1]:
        left.append(i)
    else:
        raise ValueError("unaccounted position for node %r" % (node,))

def fork_strand_info(i, node, oric, ter, left_forward, left_reverse,
        right_forward, right_reverse, rest):
    """
    Regulators are ignored for now.
    """
    if isinstance(node, Gene):
        start = node.position_start
        direction = node.strand
    elif isinstance(node, (TranscriptionUnit, Operon)):
        start = min(gene.position_start for gene in node.genes)
        direction = set(gene.strand for gene in node.genes)
        if len(direction) > 1:
            raise ValueError("conflicting strand information %r" % (node,))
        direction = direction.pop()
    else:
        rest.append(i)
        return
    if start > ter and start < oric[0]:
        if direction == "reverse":
            right_reverse.append(i)
        else:
            right_forward.append(i)
    elif start < ter or start > oric[1]:
        if direction == "reverse":
            left_reverse.append(i)
        else:
            left_forward.append(i)
    else:
        raise ValueError("unaccounted position for node %r" % (node,))

