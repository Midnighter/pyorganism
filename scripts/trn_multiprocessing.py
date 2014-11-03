#!/usr/bin/env python
# -*- coding: utf-8


from __future__ import (absolute_import, unicode_literals)

import os
import sys
import logging
import argparse
import multiprocessing
import random
from logging.config import dictConfig
from glob import glob
from random import choice

import numpy as np
import networkx as nx
import pandas as pd
from progressbar import (ProgressBar, Timer, SimpleProgress, Bar, Percentage, ETA)

import pyorganism as pyorg
import pyorganism.regulation as pyreg
from pyorganism.utils import version_from_path
from meb.utils.network.randomisation import NetworkRewiring
from meb.utils.network.subgraphs import triadic_census


LOGGER = logging.getLogger()
LOGGER.addHandler(logging.StreamHandler())


###############################################################################
# Randomisation
###############################################################################

def rewire((net, prob)):
    rnd = np.random.sample
    max_tries = 200
    nodes = sorted(net.nodes_iter())
    regulating = {node for (node, deg) in net.out_degree_iter() if deg > 0}
    regulated = set(nodes) - regulating
    edges = net.edges()
    for (u, v) in edges:
        if rnd() < prob:
            targets = list(regulated - set(net.successors_iter(u)))
            new = choice(targets)
            tries = 0
            while (net.has_edge(u, new) or net.has_edge(new, u)) and max_tries < tries:
                max_tries += 1
                new = choice(targets)
            if tries == max_tries:
                continue
            net.remove_edge(u, v)
            net.add_edge(u, new)
    return net

def switch((net, flips)):
    rewirer = NetworkRewiring()
    (rnd_net, flip_rate) = rewirer.randomise(net, flip=flips, copy=False)
    return (rnd_net, flip_rate)

def main_random(args):
    locations = sorted(glob(os.path.join(args.in_path, args.glob)))
    locations = [os.path.abspath(loc) for loc in locations]
    pool = multiprocessing.Pool(args.nproc)
    bar = ProgressBar(maxval=args.rnd_num, widgets=[Timer(), " ",
        SimpleProgress(), " ", Percentage(), " ", Bar(), " ", ETA()])
    rewire_name = "grn_rewired_{0:.1f}.pkl".format(args.prob)
    switch_name = "grn_switched_{0:d}.pkl".format(args.flip_num)
    for path in locations:
        filename = os.path.join(path, "trn.pkl")
        if not os.path.exists(filename):
            continue
        ver = version_from_path(path)
        base_path = os.path.join(args.out_path, ver)
        if not os.path.isdir(base_path):
            os.makedirs(base_path)
        LOGGER.info(ver)
        trn = pyorg.read_pickle(filename)
        # we consider null-models on the projected level since that is the level
        # on which we evaluate topological quantities
        net = pyreg.to_simple(trn.to_grn())
        if args.run_rewire:
            LOGGER.info("Rewiring with probability %.1f", args.prob)
            tasks = [(net, args.prob)] * args.rnd_num
            res_it = pool.imap_unordered(rewire, tasks)
            rands = list()
            bar.start()
            for rnd in res_it:
                rands.append(rnd)
                bar += 1
            bar.finish()
            pyorg.write_pickle(rands, os.path.join(base_path, rewire_name))
        if args.run_switch:
            LOGGER.info("Switch-randomizing each edge %d times", args.flip_num)
            tasks = [(net, args.flip_num)] * args.rnd_num
            res_it = pool.imap_unordered(switch, tasks)
            success = list()
            rands = list()
            bar.start()
            for (rnd, rate) in res_it:
                rands.append(rnd)
                success.append(rate)
                bar += 1
            bar.finish()
            pyorg.write_pickle(rands, os.path.join(base_path, switch_name))
            LOGGER.info("mean flip success rate: %.3G +/- %.3G", np.mean(success),
                    np.std(success))
    pool.close()

###############################################################################
# Analysis
###############################################################################

def stats((grn, version, description)):
    nodes = sorted(grn.nodes_iter())
    regulating = {node for (node, deg) in grn.out_degree_iter() if deg > 0}
    regulated = set(nodes) - regulating
    components = sorted(nx.weakly_connected_components(grn), key=len,
            reverse=True)
    data = dict()
    census = triadic_census(grn)
    forward = census["030T"]
    feedback = census["030C"]
    num_cycles = sum(1 for cyc in nx.simple_cycles(grn) if len(cyc) > 2)
    in_deg = [grn.in_degree(node) for node in regulated]
    out_deg = [grn.out_degree(node) for node in regulating]
    data["version"] = version
    data["num_components"] = len(components)
    data["largest_component"] = len(components[0])
    data["feed_forward"] = forward
    data["feedback"] = feedback
    data["cycles"] = num_cycles
    data["regulated_in_deg"] = np.mean(in_deg)
    data["regulating_out_deg"] = np.mean(out_deg)
    data["null_model"] = description
    stats = pd.DataFrame(data, index=[1])
    return stats

def main_analysis(args):
    locations = sorted(glob(os.path.join(args.in_path, args.glob)))
    locations = [os.path.abspath(loc) for loc in locations]
    tasks = list()
    if args.run_rewire:
        tasks.extend([(loc, "rewired", args.choose, args.prob) for loc in locations])
    if args.run_switch:
        tasks.extend([(loc, "switch", args.choose) for loc in locations])
    out_name = os.path.join(args.out_path, args.file)
    if not os.path.isdir(args.out_path):
        os.makedirs(args.out_path)
    if os.path.exists(out_name):
        result = pd.read_csv(out_name, sep=str(";"), dtype={"version": str},
                encoding=args.encoding)
    else:
        result = pd.DataFrame(columns=["version", "num_components",
            "largest_component", "feed_forward", "feedback", "cycles",
            "regulated_in_deg", "regulating_out_deg", "null_model"])
    pool = multiprocessing.Pool(args.nproc)
    rewire_name = "grn_rewired_{0:.1f}.pkl".format(args.prob)
    rewire_descr = "rewired {0:.1f}".format(args.prob)
    switch_name = "grn_switched_{0:d}.pkl".format(args.flip_num)
    switch_descr = "switch {0:d}".format(args.flip_num)
    for path in locations:
        ver = version_from_path(path)
        LOGGER.info(ver)
        if args.run_rewire:
            filename = os.path.join(path, rewire_name)
            if os.path.exists(filename):
                LOGGER.info("Analyzing rewired networks (probability %.1f)",
                        args.prob)
                nets = pyorg.read_pickle(filename)
                if args.choose is not None:
                    chosen_num = min(args.choose, len(nets))
                    LOGGER.info("Using %d/%d random networks", chosen_num,
                            len(nets))
                    nets = random.sample(nets, chosen_num)
                tasks = [(net, ver, rewire_descr) for net in nets]
                res_it = pool.imap_unordered(stats, tasks)
                frames = list()
                bar = ProgressBar(maxval=len(tasks), widgets=[Timer(), " ",
                        SimpleProgress(), " ", Percentage(), " ", Bar(), " ",
                        ETA()]).start()
                for df in res_it:
                    frames.append(df)
                    bar += 1
                bar.finish()
                result = result.append(pd.concat(frames, ignore_index=True),
                        ignore_index=True)
                result.to_csv(out_name, header=True, index=False,
                        sep=str(";"), encoding=args.encoding)
        if args.run_switch:
            filename = os.path.join(path, switch_name)
            if os.path.exists(filename):
                LOGGER.info("Analyzing switched networks (flip number %d)",
                        args.flip_num)
                nets = pyorg.read_pickle(filename)
                if args.choose is not None:
                    chosen_num = min(args.choose, len(nets))
                    LOGGER.info("Using %d/%d random networks", chosen_num, len(nets))
                    nets = random.sample(nets, chosen_num)
                tasks = [(net, ver, switch_descr) for net in nets]
                res_it = pool.imap_unordered(stats, tasks)
                frames = list()
                bar = ProgressBar(maxval=len(tasks), widgets=[Timer(), " ",
                        SimpleProgress(), " ", Percentage(), " ", Bar(), " ",
                        ETA()]).start()
                for df in res_it:
                    frames.append(df)
                    bar += 1
                bar.finish()
                result = result.append(pd.concat(frames, ignore_index=True),
                        ignore_index=True)
                result.to_csv(out_name, header=True, index=False,
                        sep=str(";"), encoding=args.encoding)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=None)
    parser.add_argument("-v", "--version", action="version", version="0.1")
    parser.add_argument("--log-level", dest="log_level", default="INFO",
            help="Log level, i.e., DEBUG, INFO, WARN, ERROR, CRITICAL (default: %(default)s)")
    parser.add_argument("--encoding", dest="encoding", default="utf-8",
            help="File encoding to assume (default: %(default)s)")
    parser.add_argument("--no-rewire", dest="run_rewire", action="store_false",
            default=True, help="Avoid creation or analysis of rewired TRNs")
    parser.add_argument("--no-switch", dest="run_switch", action="store_false",
            default=True, help="Avoid creation or analysis of switch-randomized TRNs")
    parser.add_argument("-g", "--glob", dest="glob", default="[0-9].[0-9]",
            help="Glob pattern for RegulonDB version directories (default: %(default)s)")
    parser.add_argument("-i", "--input", dest="in_path", default="RegulonDBObjects",
            help="Base directory for data input (default: %(default)s)")
    parser.add_argument("-o", "--output", dest="out_path", default="RegulonDBObjects",
            help="Base directory for data output (default: %(default)s)")
    parser.add_argument("-p", "--probability", dest="prob",
            default=0.1, type=float,
            help="Probability for rewiring a link (default: %(default)s)")
    parser.add_argument("-f", "--flip-num", dest="flip_num",
            default=int(1E02), type=int,
            help="Number of attempts to switch each link (default: %(default)s)")
    parser.add_argument("-n", "--nproc", dest="nproc",
            default=multiprocessing.cpu_count(), type=int,
            help="Number of processors to use (default: %(default)s)")
    subparsers = parser.add_subparsers(help="sub-command help")
# randomization
    parser_rnd = subparsers.add_parser("randomization",
            help="Rewire or randomize the TRN as a statistical null model")
    parser_rnd.add_argument("-r", "--rnd-num", dest="rnd_num",
            default=int(1E03), type=int,
            help="Number of rewired or randomized TRNs to generate (default: %(default)s)")
    parser_rnd.set_defaults(func=main_random)
# analysis
    parser_anal = subparsers.add_parser("analysis",
            help="Analyze the rewired or randomized TRNs")
    parser_anal.add_argument("-f", "--filename", dest="file",
            default="grn_random_stats.csv",
            help="Name of the file that statistics are written to (default: %(default)s)")
    parser_anal.add_argument("-c", "--choose", dest="choose", type=int,
            help="Size of the subset of random networks to evaluate")
    parser_anal.set_defaults(func=main_analysis)
    args = parser.parse_args()
    logger = multiprocessing.log_to_stderr()
    dictConfig({"version": 1, "incremental": True, "root": {"level":
        args.log_level}, logger.name: {"level": args.log_level}})
    try:
        sys.exit(args.func(args))
    except: # we want to catch everything
        (err, msg, trace) = sys.exc_info()
        # interrupt remote kernels and clear job queue
        raise err, msg, trace
    finally:
        logging.shutdown()

