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
        ver = os.path.basename(path)
        if not ver:
            ver = os.path.basename(os.path.dirname(path))
        base_path = os.path.join(args.out_path, ver)
        if not os.path.isdir(base_path):
            os.makedirs(base_path)
        LOGGER.info(ver)
        trn = pyorg.read_pickle(filename)
        # we consider null-models on the projected level since that is the level
        # on which we evaluate topological quantities
        net = trn.to_grn().to_simple()
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

def stats(grn, version, description):
    nodes = sorted(grn.nodes_iter())
    regulating = {node for (node, deg) in grn.out_degree_iter() if deg > 0}
    regulated = set(nodes) - regulating
    components = sorted(nx.weakly_connected_components(grn), key=len,
            reverse=True)
    data = dict()
    census = triadic_census(grn)
    forward = census["030T"]
    feedback = census["030C"]
    cycles = list(nx.simple_cycles(grn))
    in_deg = [grn.in_degree(node) for node in regulated]
    out_deg = [grn.out_degree(node) for node in regulating]
    data["version"] = version
    data["num_components"] = len(components)
    data["largest_component"] = len(components[0])
    data["feed_forward"] = forward
    data["feedback"] = feedback
    data["cycles"] = len(cycles)
    data["regulated_in_deg"] = np.mean(in_deg)
    data["regulating_out_deg"] = np.mean(out_deg)
    data["null_model"] = description
    stats = pd.DataFrame(data, index=[1])
    return stats

def err_stats(version, description):
    data = dict()
    data["version"] = version
    data["num_components"] = None
    data["largest_component"] = None
    data["feed_forward"] = None
    data["feedback"] = None
    data["cycles"] = None
    data["regulated_in_deg"] = None
    data["regulating_out_deg"] = None
    data["null_model"] = description
    return pd.DataFrame(data, index=[1])

def null_stats(parameters):
    (base_dir, task, choose) = parameters[:3]
    logger = multiprocessing.get_logger()
    ver = os.path.basename(base_dir)
    if not ver:
        ver = os.path.basename(os.path.dirname(base_dir))
    if task == "rewired":
        prob = parameters[3]
        desc = "rewired {0:.1f}".format(prob)
        filename = "trn_rewired_{0:.1f}.pkl".format(prob)
    elif task == "switch":
        desc = "switch"
        filename = "trn_random.pkl"
    try:
        nets = pyorg.read_pickle(os.path.join(base_dir, filename))
    except (OSError, IOError, EOFError):
        (err, msg, trace) = sys.exc_info()
        logger.error("Version: '%s' Task: '%s'", ver, task)
        logger.error(str(msg))
        return err_stats(ver, desc)
    chosen = random.sample(nets, choose)
    logger.info("%d/%d random networks", len(chosen), len(nets))
    return pd.concat([stats(net, ver, desc) for net in nets], ignore_index=True)

def do_task(task, results):
    loc = task[0]
    ver = os.path.basename(loc)
    if not ver:
        ver = os.path.basename(os.path.dirname(loc))
    mask = (results["version"] == ver)
    method = task[1]
    choose = task[2]
    return sum(method in descr for descr in results.loc[mask, "null_model"]) < choose

def main_analysis(args):
    locations = sorted(glob(os.path.join(args.in_path, args.glob)))
    locations = [os.path.abspath(loc) for loc in locations]
    tasks = list()
    if args.run_rewire:
        tasks.extend([(loc, "rewired", args.choose, args.prob) for loc in locations])
    if args.run_switch:
        tasks.extend([(loc, "switch", args.choose) for loc in locations])
    filename = os.path.join(args.out_path, args.file)
    if os.path.exists(filename):
        result = pd.read_csv(filename, sep=str(";"), dtype={"version": str},
                encoding=args.encoding)
        tasks = [task for task in tasks if do_task(task, result)]
    else:
        result = pd.DataFrame(columns=["version", "num_components",
            "largest_component", "feed_forward", "feedback", "cycles",
            "regulated_in_deg", "regulating_out_deg", "null_model"])
    pool = multiprocessing.Pool(args.nproc)
    res_it = pool.imap_unordered(null_stats, tasks)
    pool.close()
    bar = ProgressBar(maxval=len(tasks), widgets=[Timer(), " ",
        SimpleProgress(), " ", Percentage(), " ", Bar(), " ", ETA()]).start()
    for df in res_it:
        result = result.append(df, ignore_index=True)
        result.to_csv(filename, header=True, index=False,
                sep=str(";"), encoding=args.encoding)
        bar += 1
    bar.finish()


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
    parser_rnd.add_argument("-p", "--probability", dest="prob",
            default=0.1, type=float,
            help="Probability for rewiring a link (default: %(default)s)")
    parser_rnd.add_argument("-f", "--flip-num", dest="flip_num",
            default=int(1E02), type=int,
            help="Number of attempts to switch each link (default: %(default)s)")
    parser_rnd.set_defaults(func=main_random)
# analysis
    parser_anal = subparsers.add_parser("analysis",
            help="Analyze the rewired or randomized TRNs")
    parser_anal.add_argument("-p", "--probability", dest="prob",
            default=0.1, type=float,
            help="Probability for rewiring a link to analyze (default: %(default)s)")
    parser_anal.add_argument("-f", "--filename", dest="file",
            default="trn_random_stats.csv",
            help="Name of the file that statistics are written to (default: %(default)s)")
    parser_anal.add_argument("-c", "--choose", dest="choose",
            default=int(1E02), type=int,
            help="Size of the subset of random networks to evaluate (default: %(default)s)")
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

