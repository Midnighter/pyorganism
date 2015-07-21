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

import numpy as np
from builtins import dict

from . import shuffling
from . import continuous_wrapper as con
from .elements import (Gene, Regulator, TranscriptionUnit, Operon)
from .. import miscellaneous as misc


__all__ = [
    "ContinuousControl",
    "form_series"
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

    def from_trn(self, trn, nodes=None, node2id=None, function="regulatory"):
        """
        Prepare data structures for continuous control analysis using a
        transcriptional regulatory network.

        Parameters
        ----------
        trn: nx.(Multi-)DiGraph
            The transcriptional regulatory network (TRN).
        node2id: dict (optional)
            A mapping from nodes in the network to indices running from 0 to (N - 1).
        function: hashable (optional)
            In case of a DiGraph, the keyword for edge data that describes the
            regulatory function. For a MultiDiGraph the key is used.
        """
        if len(trn) < 2 or trn.size() < 1:
            raise ValueError("aborting due to small network size")
        if nodes is None:
            self.nodes = sorted(trn.nodes())
        else:
            self.nodes = nodes
        self.num_nodes = len(self.nodes)
        if node2id is None:
            self.node2id = {n: i for (i, n) in enumerate(self.nodes)}
        else:
            self.node2id = node2id
        sources = list()
        targets = list()
        functions = list()
        if trn.is_multigraph():
            edge_iter = ((u, v, k) for (u, v, k) in trn.edges_iter(keys=True))
        else:
            edge_iter = ((u, v, data[function]) for (u, v, data) in\
                    trn.edges_iter(data=True))
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

    def from_gpn(self, gpn, nodes=None, node2id=None):
        """
        Prepare data structures for continuous control analysis using a
        gene proximity network.

        Parameters
        ----------
        gpn: nx.Graph
            The gene proximity network (GPN).
        node2id: dict (optional)
            A mapping from nodes in the network to indices running from 0 to (N - 1).
        """
        if len(gpn) < 2 or gpn.size() < 1:
            raise ValueError("aborting due to small network size")
        if nodes is None:
            self.nodes = sorted(gpn.nodes())
        else:
            self.nodes = nodes
        self.num_nodes = len(self.nodes)
        if node2id is None:
            self.node2id = {n: i for (i, n) in enumerate(self.nodes)}
        else:
            self.node2id = node2id
        self.num_links = gpn.size()
        self.sources = np.zeros(self.num_links, dtype=np.int32)
        self.targets = np.zeros(self.num_links, dtype=np.int32)
        self.functions = None
        for (i, (u, v)) in enumerate(gpn.edges_iter()):
            self.sources[i] = self.node2id[u]
            self.targets[i] = self.node2id[v]

    def series_ctc(self, expression, sampling, measure, random_num=1E04, delay=0,
            **extra_args):
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
            sample = self._sampling[sampling]
        except KeyError:
            raise ValueError("'{0}' is an unknown sampling"\
                    " method.".format(sampling))
        try:
            func = self._measures[measure]
        except KeyError:
            raise ValueError("'{0}' is an unknown control measure."\
                    " Please try one of:\n{1}".format(measure,
                    "\n* ".join(self._measures.keys())))
        if "comparison" in measure:
            expression = np.diff(expression, axis=0) # expression was transposed for C code
        num_points = expression.shape[0]
        num_features = expression.shape[1]
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
        for (i, rnd_series) in enumerate(sample(expression, random_num,
            **extra_args)):
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

def form_series(net, df, feature2node, node2feature=None, nodes=None,
        include_fork=False, include_strand=False, oric=(3923767, 3923998),
        size=4639675):
    """
    Convert a pandas DataFrame to an expression matrix.
    """
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
    sub = net.subgraph(active)
    nom = len(active)
    denom = len(nodes)
    LOGGER.info("%d/%d active nodes (%.2f%%)", nom, denom, nom / denom * 100.0)
    expression = np.zeros((df.shape[1], len(active)))
    if include_fork:
        ter = int((np.mean(oric) + size / 2) % size)
        rest = list()
        if include_strand:
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
            extra = dict(blocks=np.cumsum(blocks, dtype=np.int32))
            # organize expression into complete blocks
            expression = np.ascontiguousarray(expression[:, left_forward + left_reverse +
                right_forward + right_reverse + rest])
        else:
            left = list()
            right = list()
            for (i, node) in enumerate(active):
                expression[:, i] = data[node]
                fork_info(i, node, oric, ter, left, right, rest)
            blocks = [0, len(left), len(right), len(rest)]
            extra = dict(blocks=np.cumsum(blocks, dtype=np.int32))
            # organize expression into complete blocks
            expression = np.ascontiguousarray(expression[:, left + right + rest])
    else:
        for (i, node) in enumerate(active):
            expression[:, i] = data[node]
        extra = dict()
    expression[~np.isfinite(expression)] = 0.0
    return (sub, expression, df.columns, extra)

