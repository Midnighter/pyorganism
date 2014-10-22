#!/usr/bin/env python
# -*- coding: utf-8


from __future__ import (absolute_import, unicode_literals)

import os
import sys
import logging
import argparse
from logging.config import dictConfig
from glob import glob
from random import choice

import numpy as np
import networkx as nx
import pandas as pd
from IPython.parallel import (interactive, Client)
from progressbar import (ProgressBar, Timer, SimpleProgress, Bar, Percentage, ETA)

import pyorganism as pyorg
from pyorganism.regulation import trn2grn
from meb.utils.network.randomisation import NetworkRewiring
from meb.utils.network.subgraphs import triadic_census


LOGGER = logging.getLogger()
LOGGER.addHandler(logging.StreamHandler())


###############################################################################
# Randomisation
###############################################################################

def load_data(locations):
    tr_nets = dict()
    versions = list()
    drop = list()
    for path in locations:
        ver = os.path.basename(path)
        if not ver:
            ver = os.path.basename(os.path.dirname(path))
        try:
            pyorg.read_pickle(os.path.join(path, "genes.pkl"))
            pyorg.read_pickle(os.path.join(path, "transcription_factors.pkl"))
            trn = pyorg.read_pickle(os.path.join(path, "trn.pkl"))
            tr_nets[ver] = trn
            versions.append(ver)
        except IOError:
            drop.append(path)
            continue
    for path in drop:
        locations.remove(path)
    return (tr_nets, versions)

@interactive
def rewire(version):
    net = globals()["TRN"][version].copy()
    prob = globals()["prob"]
    rnd = np.random.sample
    nodes = sorted(net.nodes_iter())
    regulating = {node for (node, deg) in net.out_degree_iter() if deg > 0}
    regulated = set(nodes) - regulating
    edges = net.edges(data=True, keys=True)
    for (u, v, key, data) in edges:
        if rnd() < prob:
            targets = list(regulated - set(net.successors(u)))
            new = choice(targets)
            while net.has_edge(u, new, key):
                new = choice(targets)
            net.remove_edge(u, v, key)
            net.add_edge(u, new, key=key, **data)
    return net

def rewiring(lb_view, versions, args):
    bar = ProgressBar(maxval=args.rnd_num, widgets=[Timer(), " ",
        SimpleProgress(), " ", Percentage(), " ", Bar(), " ", ETA()])
    rands = list()
    for ver in versions:
        LOGGER.info(ver)
        res_it = lb_view.map(rewire, [ver] * args.rnd_num, block=False, ordered=False)
        bar.start()
        for rng in res_it:
            rands.append(rng)
            bar += 1
        bar.finish()
        pyorg.write_pickle(rands, os.path.join(args.out_path, ver,
                "trn_rewired_{0:.1f}.pkl".format(args.prob)))
        lb_view.purge_results("all")
        del rands[:] # activate garbage collection in loop

@interactive
def null_model(version):
    net = nx.DiGraph(globals()["TRN"][version])
    flips = globals()["flip_num"]
    nx.convert_node_labels_to_integers(net, ordering="sorted",
            label_attribute="element")
    rewirer = NetworkRewiring()
    (rnd_net, flip_rate) = rewirer.randomise(net, flip=flips, copy=False)
    return (rnd_net, flip_rate)

def randomisation(lb_view, versions, args):
    bar = ProgressBar(maxval=args.rnd_num, widgets=[Timer(), " ",
        SimpleProgress(), " ", Percentage(), " ", Bar(), " ", ETA()])
    for ver in versions:
        LOGGER.info(ver)
        res_it = lb_view.map(null_model, [ver] * args.rnd_num, block=False, ordered=False)
        bar.start()
        rands = list()
        success = list()
        for (rnd_net, flip_rate) in res_it:
            rands.append(rnd_net)
            success.append(flip_rate)
            bar +=1
        bar.finish()
        lb_view.purge_results("all")
        LOGGER.info("mean flip success rate: %.3G +/- %.3G", np.mean(success),
                np.std(success))
        pyorg.write_pickle(rands, os.path.join(args.out_path, ver,
                "trn_random.pkl"))
        del rands[:] # activate garbage collection in loop

def main_random(rc, args):
    locations = sorted(glob(os.path.join(args.in_path, args.glob)))
    locations = [os.path.abspath(loc) for loc in locations]
    LOGGER.info("loading data")
    (tr_nets, versions) = load_data(locations)
    LOGGER.info("remote preparation")
    dv = rc.direct_view()
    dv.execute("import os;"\
            "from random import choice;"\
            "import numpy as np;"\
            "import networkx as nx;"\
            "import pyorganism as pyorg;"\
            "from meb.utils.network.randomisation import NetworkRewiring"\
            "import logging; from IPython.config import Application;"\
            "LOGGER = Application.instance().log;"\
            "LOGGER.setLevel(logging.{level});".format(level=args.log_level),
            block=True)
    dv.push({"load_data": load_data, "locations": locations}, block=True)
    dv.execute("(TRN, versions) = load_data(locations);", block=True)
    lv = rc.load_balanced_view()
    if args.run_rewire:
        LOGGER.info("rewiring")
        dv.push({"rewire": rewire, "prob": args.prob}, block=True)
        rewiring(lv, versions, args)
    if args.run_rnd:
        LOGGER.info("randomisation")
        dv.push({"null_model": null_model, "flip_num": args.flip_num}, block=True)
        randomisation(lv, versions, args)

###############################################################################
# Analysis
###############################################################################

@interactive
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

@interactive
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

@interactive
def null_stats(base_dir, task):
    glbls = globals()
    prob = glbls["prob"]
    choose = glbls["choose"]
    logger = glbls["LOGGER"]
    ver = os.path.basename(base_dir)
    if not ver:
        ver = os.path.basename(os.path.dirname(base_dir))
    if task == "rewired":
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
    nets = [trn2grn(net) for net in chosen]
    return pd.concat([stats(net, ver, desc) for net in nets], ignore_index=True)

def filter_tasks(locations, tasks, results, num_expect=1000):
    if results.empty:
        return ([loc for loc in locations for method in tasks],
                [method for loc in locations for method in tasks])
    paths = list()
    methods = list()
    for loc in locations:
        ver = os.path.basename(loc)
        if not ver:
            ver = os.path.basename(os.path.dirname(loc))
        mask = (results["version"] == ver)
        for t in tasks:
            if sum(t in dscr for dscr in results.loc[mask,
                "null_model"]) != num_expect:
                paths.append(loc)
                methods.append(t)
    return (paths, methods)

def main_analysis(rc, args):
    locations = sorted(glob(os.path.join(args.in_path, args.glob)))
    locations = [os.path.abspath(loc) for loc in locations]
    LOGGER.info("remote preparation")
    dv = rc.direct_view()
    dv.execute("import os;"\
            "import sys;"\
            "import random;"\
            "import numpy as np;"\
            "import networkx as nx;"\
            "import pandas as pd;"\
            "from pyorganism.regulation import trn2grn;"\
            "from meb.utils.network.subgraphs import triadic_census;"\
            "import pyorganism as pyorg;"\
            "import logging; from IPython.config import Application;"\
            "LOGGER = Application.instance().log;"\
            "LOGGER.setLevel(logging.{level});".format(level=args.log_level),
            block=True)
    dv.push({"stats": stats, "err_stats": err_stats, "choose": args.choose}, block=True)
    tasks = list()
    if args.run_rewire:
        dv.push({"prob": args.prob}, block=True)
        tasks.append("rewired")
    if args.run_rnd:
        tasks.append("switch")
    filename = os.path.join(args.out_path, args.file)
    if os.path.exists(filename):
        result = pd.read_csv(filename, sep=str(";"), dtype={"version": str},
                encoding=args.encoding)
    else:
        result = pd.DataFrame(columns=["version", "num_components",
            "largest_component", "feed_forward", "feedback", "cycles",
            "regulated_in_deg", "regulating_out_deg", "null_model"])
    (locations, methods) = filter_tasks(locations, tasks, result)
    lv = rc.load_balanced_view()
    res_it = lv.map(null_stats, locations, methods, block=False, ordered=False)
    bar = ProgressBar(maxval=len(locations) * 2, widgets=[Timer(), " ",
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
    parser.add_argument("--profile", dest="profile", default="default",
            help="IPython profile to connect to cluster (default: %(default)s)")
    parser.add_argument("--cluster-id", dest="cluster_id", default=None,
            help="IPython cluster-id to connect to (default: %(default)s)")
    parser.add_argument("--log-level", dest="log_level", default="INFO",
            help="Log level, i.e., DEBUG, INFO, WARN, ERROR, CRITICAL (default: %(default)s)")
    parser.add_argument("--encoding", dest="encoding", default="utf-8",
            help="File encoding to assume (default: %(default)s)")
    parser.add_argument("--no-rewire", dest="run_rewire", action="store_false",
            default=True, help="Avoid creation or analysis of rewired TRNs")
    parser.add_argument("--no-randomize", dest="run_rnd", action="store_false",
            default=True, help="Avoid creation or analysis of randomized TRNs")
    parser.add_argument("-g", "--glob", dest="glob", default="[0-9].[0-9]",
            help="Glob pattern for RegulonDB version directories (default: %(default)s)")
    parser.add_argument("-i", "--input", dest="in_path", default="RegulonDBObjects",
            help="Base directory for data input (default: %(default)s)")
    parser.add_argument("-o", "--output", dest="out_path", default="RegulonDBObjects",
            help="Base directory for data output (default: %(default)s)")
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
    dictConfig({"version": 1, "incremental": True, "root": {"level": args.log_level}})
    remote_client = Client(profile=args.profile, cluster_id=args.cluster_id)
    try:
        sys.exit(args.func(remote_client, args))
    except: # we want to catch everything
        (err, msg, trace) = sys.exc_info()
        # interrupt remote kernels and clear job queue
        raise err, msg, trace
    finally:
        logging.shutdown()

