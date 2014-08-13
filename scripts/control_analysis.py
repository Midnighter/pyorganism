#!/usr/bin/env python
# -*- coding: utf-8 -*-


from __future__ import (division, print_function)


import sys
import os
import logging
import argparse
import json
import codecs
import shelve
import pickle

import numpy as np

import pyorganism as pyorg
import pyorganism.regulation as pyreg
import pyorganism.io.microarray as pymicro

from itertools import (izip, chain)
from signal import SIGINT
from exceptions import SystemExit
from logging.config import dictConfig

from IPython.parallel import (Client, interactive, RemoteError)
from progressbar import (ProgressBar, Timer, Bar, Percentage, ETA)

from pyorganism.io.hdf5 import ResultManager


LOGGER = logging.getLogger()
LOGGER.addHandler(logging.StreamHandler())


def check_path(path):
    if not os.path.isfile(path):
        raise IOError("file does not exist '{path}'".format(path=path))

def string_pair(key, value, level, sep=": "):
    return "{0}{1}{2}".format(str(key), sep, str(value)).rjust(level)

def print_dict(dic, level=0):
    message = list()
    for (key, value) in dic.iteritems():
        if isinstance(value, dict):
            message.append(string_pair(key, "{", level))
            message.extend(print_dict(value, level + 2))
            message.append(string_pair("", "}", level, sep=""))
        elif isinstance(value, list):
            message.append(string_pair(key, "[...]", level))
        elif isinstance(value, set):
            message.append(string_pair(key, "{...}", level))
        else:
            message.append(string_pair(key, value, level))
    return message


##############################################################################
# Discrete
##############################################################################


def load_discrete(organism, config):
    LOGGER.info("Loading differentially expressed genes:")
    experiments = config["experiments"]
    for (filename, name, reader, extra) in izip(experiments["paths"],
            experiments["names"], experiments["readers"],
            experiments["args"]):
        reader_func = getattr(pymicro, reader)
        path = os.path.join(experiments["base"], filename)
        check_path(path)
        LOGGER.info("  %s: '%s'", name, path)
        organism.activity[name] = reader_func(path)

def simple_discrete(control_type, df, name2gene):
    eligible = df["name"].notnull()
    active = [name2gene[name] for name in df["name"][eligible]\
            if name2gene[name] is not None]
    LOGGER.info("  mapped %d/%d active genes (%.2f%%)", len(active), len(df),
            len(active) / len(df) * 100.0)
    results = dict()
    results["gene"] = active
    if control_type == "digital":
        results["gene"] = pyreg.active_genes_and_tf(active)
    results["tu"] = pyreg.active_tu(active)
    results["operon"] = pyreg.active_operons(active)
    return results

def ratio_discrete(control_type, df, name2gene):
    eligible = df["name"].notnull()
    up = (df["ratio (A/B)"] > 1.0) & eligible
    down = (df["ratio (A/B)"] < 1.0) & eligible
    up_reg = [name2gene[name] for name in df["name"][up]\
            if name2gene[name] is not None]
    total = up.sum()
    LOGGER.info("  mapped %d/%d up-regulated genes (%.2f%%)", len(up_reg), total,
            len(up_reg) / total * 100.0)
    down_reg = [name2gene[name] for name in df["name"][down]\
            if name2gene[name] is not None]
    total = down.sum()
    LOGGER.info("  mapped %d/%d down-regulated genes (%.2f%%)", len(down_reg),
            total, len(down_reg) / total * 100.0)
    results = dict()
    results["gene"] = dict()
    results["gene"]["up"] = up_reg
    results["gene"]["down"] = down_reg
    results["tu"] = dict()
    results["tu"]["up"] = pyreg.active_tu(up_reg)
    results["tu"]["down"] = pyreg.active_tu(down_reg)
    results["operon"] = dict()
    results["operon"]["up"] = pyreg.active_operons(up_reg)
    results["operon"]["down"] = pyreg.active_operons(down_reg)
    return results

def discrete_jobs(organism, config, *args):
    LOGGER.info("Generating discrete job specifications:")
    jobs = list()
    analysis = config["analysis"]
    for version in config["versions"]:
        for (cntrl_name, experiments, setups,
                control, ctc, measures, random_num, robustness_num,
                rob_extra, projections) in izip(analysis["control_types"],
                analysis["experimental_sets"], analysis["experimental_setups"],
                analysis["control"], analysis["ctc"],
                analysis["measures"], analysis["random_num"],
                analysis["robustness_num"], analysis["robustness_args"],
                config["network"]["projections"]):
            for method in ctc:
                for basis in projections:
                    for ms_name in measures:
                        for (exp_name, exp_setup) in izip(experiments, setups):
                            if exp_setup == "ratio_discrete":
                                for direction in ["up", "down"]:
                                    spec = dict()
                                    spec["version"] = version
                                    spec["continuous"] = config["continuous"]
                                    spec["control_type"] = cntrl_name
                                    spec["experiment"] = exp_name
                                    spec["projection"] = basis
                                    spec["setup"] = exp_setup
                                    spec["direction"] = direction
                                    spec["control"] = control
                                    spec["ctc"] = method
                                    spec["measure"] = ms_name
                                    spec["random_num"] = random_num
                                    spec["robustness_num"] = robustness_num
                                    spec["robustness_args"] = rob_extra
                                    jobs.append(spec)
                            else:
                                spec = dict()
                                spec["version"] = version
                                spec["continuous"] = config["continuous"]
                                spec["control_type"] = cntrl_name
                                spec["experiment"] = exp_name
                                spec["projection"] = basis
                                spec["setup"] = exp_setup
                                spec["control"] = control
                                spec["ctc"] = method
                                spec["measure"] = ms_name
                                spec["random_num"] = random_num
                                spec["robustness_num"] = robustness_num
                                spec["robustness_args"] = rob_extra
                                jobs.append(spec)
    LOGGER.info("  %d jobs", len(jobs))
    return jobs

@interactive
def discrete_worker(spec):
    LOGGER.debug(spec)
    version = spec["version"]
    cntrl_type = spec["control_type"]
    global_vars = globals()
    control = getattr(pyreg, spec["control"])
    ctc = getattr(pyreg, spec["ctc"])
    measure = getattr(pyreg, spec["measure"])
    net = global_vars["networks"][version][cntrl_type][spec["projection"]]
    if cntrl_type == "analog":
        active = global_vars["prepared"][version][cntrl_type][spec["experiment"]][spec["projection"]][spec["direction"]]
    else:
        active = global_vars["prepared"][version][cntrl_type][spec["experiment"]][spec["projection"]]
    LOGGER.debug(len(active))
    effective = pyreg.effective_network(net, active)
    res_cntrl = control(effective, measure=measure)
    (res_ctc, samples) = ctc(effective, net, measure=measure,
            random_num=spec["random_num"], return_sample=True)
    return (spec, res_cntrl, res_ctc, samples)

def discrete_result(manager, spec, res_cntrl, res_ctc, samples):
    manager.append(version=spec["version"], control_type=spec["control_type"],
            continuous=spec["continuous"], strain=spec["experiment"],
            projection=spec["projection"], setup=spec["setup"], control_strength=res_cntrl,
            control_method=spec["control"], ctc=res_ctc, ctc_method=spec["ctc"],
            measure=spec["measure"], direction=spec.get("direction", None),
            samples=samples)


##############################################################################
# Continuous
##############################################################################


def simple_continuous(control_type, df, feature2gene):
    df[df <= 0.0] = np.nan
    df = df.apply(pyorg.norm_zero2unity, axis=1, raw=True)
    times = list()
    actives = list()
    levels = list()
    for col in df.columns:
        eligible = df[col][np.isfinite(df[col])]
        mask = [feature2gene[name] is not None for name in eligible.index]
        active = [feature2gene[name] for name in eligible[mask].index]
        level = eligible[mask]
        LOGGER.info("        %s min: %d/%d active genes (%.2f%%)", col,
                len(active), len(df), len(active) / len(df) * 100.0)
        times.append(col)
        actives.append(active)
        levels.append(level)
    results = dict()
    results["time"] = times
    results["gene"] = dict()
    results["gene"]["active"] = actives
    results["gene"]["levels"] = levels
    if control_type == "digital":
        trn_actives = list()
        trn_levels = list()
        for i in range(len(actives)):
            (active, level) = pyreg.gene_and_tf_levels(actives[i], levels[i])
            trn_actives.append(active)
            trn_levels.append(level)
        results["gene"]["active"] = trn_actives
        results["gene"]["levels"] = trn_levels
    tu_actives = list()
    tu_levels = list()
    for i in range(len(actives)):
        (active, level) = pyreg.tu_levels(actives[i], levels[i])
        tu_actives.append(active)
        tu_levels.append(level)
    results["tu"] = dict()
    results["tu"]["active"] = tu_actives
    results["tu"]["levels"] = tu_levels
    op_actives = list()
    op_levels = list()
    for i in range(len(actives)):
        (active, level) = pyreg.operon_levels(actives[i], levels[i])
        op_actives.append(active)
        op_levels.append(level)
    results["operon"] = dict()
    results["operon"]["active"] = op_actives
    results["operon"]["levels"] = op_levels
    return results

def rate_continuous(results):
    times = results["time"]
    actives = results["gene"]["active"]
    levels = results["gene"]["levels"]
    rate_time = list()
    rate_actives = list()
    rate_levels = list()
    for i in range(len(times) - 1):
        t1 = dict(izip(actives[i], levels[i]))
        t2 = dict(izip(actives[i + 1], levels[i + 1]))
        joint = sorted(set(actives[i]) & set(actives[i + 1]))
        active = list()
        level = list()
        for gene in joint:
            if gene in t1 and gene in t2:
                active.append(gene)
                level.append(t2[gene] - t1[gene])
        LOGGER.info("        %s - %s min: %d active genes", times[i],
                times[i + 1], len(active))
        rate_time.append(times[i + 1])
        rate_actives.append(active)
        rate_levels.append(level)
    results["rate"] = dict()
    results["rate"]["time"] = rate_time
    results["rate"]["gene"] = dict()
    results["rate"]["gene"]["active"] = rate_actives
    results["rate"]["gene"]["levels"] = rate_levels
    # TU
    actives = results["tu"]["active"]
    levels = results["tu"]["levels"]
    rate_time = list()
    rate_actives = list()
    rate_levels = list()
    for i in range(len(times) - 1):
        t1 = dict(izip(actives[i], levels[i]))
        t2 = dict(izip(actives[i + 1], levels[i + 1]))
        joint = sorted(set(actives[i]) & set(actives[i + 1]))
        active = list()
        level = list()
        for tu in joint:
            if tu in t1 and tu in t2:
                active.append(tu)
                level.append(t2[tu] - t1[tu])
        LOGGER.info("        %s - %s min: %d active TUs", times[i],
                times[i + 1], len(joint))
        rate_time.append(times[i + 1])
        rate_actives.append(active)
        rate_levels.append(level)
    results["rate"]["tu"] = dict()
    results["rate"]["tu"]["active"] = rate_actives
    results["rate"]["tu"]["levels"] = rate_levels
    # operon
    actives = results["operon"]["active"]
    levels = results["operon"]["levels"]
    rate_time = list()
    rate_actives = list()
    rate_levels = list()
    for i in range(len(times) - 1):
        t1 = dict(izip(actives[i], levels[i]))
        t2 = dict(izip(actives[i + 1], levels[i + 1]))
        joint = sorted(set(actives[i]) & set(actives[i + 1]))
        active = list()
        level = list()
        for op in joint:
            if op in t1 and op in t2:
                active.append(op)
                level.append(t2[op] - t1[op])
        LOGGER.info("        %s - %s min: %d active operons", times[i],
                times[i + 1], len(active))
        rate_time.append(times[i + 1])
        rate_actives.append(active)
        rate_levels.append(level)
    results["rate"]["operon"] = dict()
    results["rate"]["operon"]["active"] = rate_actives
    results["rate"]["operon"]["levels"] = rate_levels

def delayed_continuous(results, delta):
    times = results["time"]
    actives = results["gene"]["active"]
    levels = results["gene"]["levels"]
    delayed_actives = list()
    delayed_levels = list()
    delta_levels = list()
    for i in range(len(times)):
        j = (i + delta) % len(times)
        t1 = dict(izip(actives[i], levels[i]))
        t2 = dict(izip(actives[j], levels[j]))
        joint = sorted(set(actives[i]) & set(actives[j]))
        active = list()
        level = list()
        delayed = list()
        for gene in joint:
            if gene in t1 and gene in t2:
                active.append(gene)
                level.append(t1[gene])
                delayed.append(t2[gene])
        LOGGER.info("        %s -> %s min: %d active genes", times[i],
                times[j], len(joint))
        delayed_actives.append(active)
        delayed_levels.append(level)
        delta_levels.append(delayed)
    results["delayed"][delta] = dict()
    results["delayed"][delta]["gene"] = dict()
    results["delayed"][delta]["gene"]["active"] = delayed_actives
    results["delayed"][delta]["gene"]["levels"] = delayed_levels
    results["delayed"][delta]["gene"]["delayed"] = delta_levels
    # TU
    actives = results["tu"]["active"]
    levels = results["tu"]["levels"]
    delayed_actives = list()
    delayed_levels = list()
    delta_levels = list()
    for i in range(len(times)):
        j = (i + delta) % len(times)
        t1 = dict(izip(actives[i], levels[i]))
        t2 = dict(izip(actives[j], levels[j]))
        joint = sorted(set(actives[i]) & set(actives[j]))
        active = list()
        level = list()
        delayed = list()
        for tu in joint:
            if tu in t1 and tu in t2:
                active.append(tu)
                level.append(t1[tu])
                delayed.append(t2[tu])
        LOGGER.info("        %s -> %s min: %d active TUs", times[i],
                times[j], len(joint))
        delayed_actives.append(active)
        delayed_levels.append(level)
        delta_levels.append(delayed)
    results["delayed"][delta]["tu"] = dict()
    results["delayed"][delta]["tu"]["active"] = delayed_actives
    results["delayed"][delta]["tu"]["levels"] = delayed_levels
    results["delayed"][delta]["tu"]["delayed"] = delta_levels
    # operon
    actives = results["operon"]["active"]
    levels = results["operon"]["levels"]
    delayed_actives = list()
    delayed_levels = list()
    delta_levels = list()
    for i in range(len(times)):
        j = (i + delta) % len(times)
        t1 = dict(izip(actives[i], levels[i]))
        t2 = dict(izip(actives[j], levels[j]))
        joint = sorted(set(actives[i]) & set(actives[j]))
        active = list()
        level = list()
        delayed = list()
        for op in joint:
            if op in t1 and op in t2:
                active.append(op)
                level.append(t1[op])
                delayed.append(t2[op])
        LOGGER.info("        %s -> %s min: %d active operons", times[i],
                times[j], len(joint))
        delayed_actives.append(active)
        delayed_levels.append(level)
        delta_levels.append(delayed)
    results["delayed"][delta]["operon"] = dict()
    results["delayed"][delta]["operon"]["active"] = delayed_actives
    results["delayed"][delta]["operon"]["levels"] = delayed_levels
    results["delayed"][delta]["operon"]["delayed"] = delta_levels

def shuffle(df, n_times=1, axis=0):
    df = df.copy()
    axis ^= 1 # make axes between numpy and pandas behave as expected (only 2D!)
    for _ in range(n_times):
        for view in np.rollaxis(df.values, axis):
            np.random.shuffle(view)
    return df

def randomised_continuous(control_type, df, feature2gene):
    df = shuffle(df)
    return simple_continuous(control_type, df, feature2gene)

def fully_randomised_continuous(control_type, df, feature2gene):
    df = shuffle(df)
    df = shuffle(df, axis=1)
    return simple_continuous(df, feature2gene)

def load_continuous(organism, config):
    LOGGER.info("Loading continuous expression data:")
    experiments = config["experiments"]
    for (filename, name, reader, extra) in izip(experiments["paths"],
            experiments["names"], experiments["readers"], experiments["args"]):
        reader_func = getattr(pymicro, reader)
        path = os.path.join(experiments["base"], filename)
        check_path(path)
        LOGGER.info("  %s: '%s'", name, path)
        organism.activity[name] = reader_func(path, mutants=extra)

def continuous_jobs(organism, config, namespace):
    LOGGER.info("Generating continuous job specifications:")
    jobs = list()
    analysis = config["analysis"]
    for version in config["versions"]:
        for (cntrl_type, experiments, setups, control, ctc, measures, random_num,
                robustness_num, rob_extra, projections) in izip(analysis["control_types"],
                analysis["experimental_sets"], analysis["experimental_setups"],
                analysis["control"], analysis["ctc"], analysis["measures"],
                analysis["random_num"], analysis["robustness_num"],
                analysis["robustness_args"], config["network"]["projections"]):
            for method in ctc:
                for basis in projections:
                    for ms_name in measures:
                        if method.startswith("delayed"):
                            for delta in analysis["delays"]:
                                for (exp_name, exp_setup) in izip(experiments, setups):
                                    times = namespace["prepared"][version][cntrl_type][exp_name]["time"]
                                    for time_point in times:
                                        spec = dict()
                                        spec["version"] = version
                                        spec["continuous"] = config["continuous"]
                                        spec["control_type"] = cntrl_type
                                        spec["time"] = time_point
                                        spec["experiment"] = exp_name
                                        spec["projection"] = basis
                                        spec["setup"] = exp_setup
                                        spec["control"] = control
                                        spec["ctc"] = method
                                        spec["measure"] = ms_name
                                        spec["delay"] = delta
                                        spec["random_num"] = random_num
                                        spec["robustness_num"] = robustness_num
                                        spec["robustness_args"] = rob_extra
                                        jobs.append(spec)
                        elif ms_name.endswith("comparison"):
                            for (exp_name, exp_setup) in izip(experiments, setups):
                                times = namespace["prepared"][version][cntrl_type][exp_name]["rate"]["time"]
                                for time_point in times:
                                    spec = dict()
                                    spec["version"] = version
                                    spec["continuous"] = config["continuous"]
                                    spec["control_type"] = cntrl_type
                                    spec["time"] = time_point
                                    spec["experiment"] = exp_name
                                    spec["projection"] = basis
                                    spec["setup"] = exp_setup
                                    spec["control"] = control
                                    spec["ctc"] = method
                                    spec["measure"] = ms_name
                                    spec["random_num"] = random_num
                                    spec["robustness_num"] = robustness_num
                                    spec["robustness_args"] = rob_extra
                                    jobs.append(spec)
                        else:
                            for (exp_name, exp_setup) in izip(experiments, setups):
                                times = namespace["prepared"][version][cntrl_type][exp_name]["time"]
                                for time_point in times:
                                    spec = dict()
                                    spec["version"] = version
                                    spec["continuous"] = config["continuous"]
                                    spec["control_type"] = cntrl_type
                                    spec["time"] = time_point
                                    spec["experiment"] = exp_name
                                    spec["projection"] = basis
                                    spec["setup"] = exp_setup
                                    spec["control"] = control
                                    spec["ctc"] = method
                                    spec["measure"] = ms_name
                                    spec["random_num"] = random_num
                                    spec["robustness_num"] = robustness_num
                                    spec["robustness_args"] = rob_extra
                                    jobs.append(spec)
    LOGGER.info("  %d jobs", len(jobs))
    return jobs

@interactive
def continuous_worker(spec):
    LOGGER.debug(spec)
    version = spec["version"]
    cntrl_type = spec["control_type"]
    global_vars = globals()
    control = getattr(pyreg, spec["control"])
    ctc = getattr(pyreg, spec["ctc"])
    measure = getattr(pyreg, spec["measure"])
    net = global_vars["networks"][version][cntrl_type][spec["projection"]]
    prepared = global_vars["prepared"][version][cntrl_type][spec["experiment"]]
    if spec["ctc"].startswith("delayed"):
        index = prepared["time"].index(spec["time"])
        active = prepared["delayed"][spec["delay"]][spec["projection"]]["active"][index]
        levels = prepared["delayed"][spec["delay"]][spec["projection"]]["levels"][index]
        delayed = prepared["delayed"][spec["delay"]][spec["projection"]]["delayed"][index]
        effective = pyreg.effective_network(net, active)
        res_cntrl = control(effective, active, levels, delayed, measure=measure)
        (res_ctc, samples) = ctc(effective, active, levels, delayed,
                random_num=spec["random_num"], measure=measure, return_sample=True)
    elif spec["measure"].endswith("comparison"):
        index = prepared["rate"]["time"].index(spec["time"])
        active = prepared["rate"][spec["projection"]]["active"][index]
        levels = prepared["rate"][spec["projection"]]["levels"][index]
        effective = pyreg.effective_network(net, active)
        res_cntrl = control(effective, active, levels, measure=measure)
        (res_ctc, samples) = ctc(effective, active, levels,
                random_num=spec["random_num"], measure=measure, return_sample=True)
    else:
        index = prepared["time"].index(spec["time"])
        active = prepared[spec["projection"]]["active"][index]
        levels = prepared[spec["projection"]]["levels"][index]
        LOGGER.debug(len(active))
        effective = pyreg.effective_network(net, active)
        res_cntrl = control(effective, active, levels, measure=measure)
        (res_ctc, samples) = ctc(effective, active, levels, random_num=spec["random_num"],
                measure=measure, return_sample=True)
    return (spec, res_cntrl, res_ctc, samples)

def continuous_result(manager, spec, res_cntrl, res_ctc, samples):
    manager.append(version=spec["version"], control_type=spec["control_type"],
            continuous=spec["continuous"], strain=spec["experiment"],
            projection=spec["projection"], setup=spec["setup"],
            control_strength=res_cntrl, control_method=spec["control"],
            ctc=res_ctc, ctc_method=spec["ctc"],
            measure=spec["measure"], time=int(float(spec["time"])),
            samples=samples, delay=spec.get("delay"))


##############################################################################
# Main
##############################################################################


def main(remote_client, args):
    config = json.load(codecs.open(args.config, encoding=args.encoding, mode="rb"))
    if config["continuous"]:
        load_func = load_continuous
        table_key = "/Continuous"
        job_gen = continuous_jobs
        worker = continuous_worker
        result = continuous_result
    else:
        load_func = load_discrete
        table_key = "/Discrete"
        job_gen = discrete_jobs
        worker = discrete_worker
        result = discrete_result
    organism = pyorg.Organism(name=config["organism"])
    load_func(organism, config)
    LOGGER.info("Load data")
    glob_vars = globals()
    data = config["data"]
    network = config["network"]
    analysis = config["analysis"]
    namespace = dict()
    namespace["genes"] = dict()
    namespace["networks"] = dict()
    namespace["prepared"] = dict()
    for version in config["versions"]:
        LOGGER.info("{0:*^78s}".format(version))
        namespace["genes"][version] = pyorg.read_pickle(os.path.join(
                data["base"], version, data["gene_path"]))
        id2gene = pyorg.read_pickle(os.path.join(data["base"], version,
                data["mapping_path"]))
        namespace["networks"][version] = dict()
        for (cntrl_type, net_file, projections) in izip(analysis["control_types"],
                network["paths"], network["projections"]):
            net = pyorg.read_pickle(os.path.join(data["base"], version, net_file))
            namespace["networks"][version][cntrl_type] = dict()
            for basis in projections:
                if basis == "gene":
                    namespace["networks"][version][cntrl_type][basis] = net
                elif basis == "tu":
                    namespace["networks"][version][cntrl_type][basis] =\
                            pyreg.to_transcription_unit_based(net)
                elif basis == "operon":
                    namespace["networks"][version][cntrl_type][basis] =\
                            pyreg.to_operon_based(net)
        namespace["prepared"][version] = dict()
        for (cntrl_type, experiments, setups) in izip(analysis["control_types"],
                 analysis["experimental_sets"], analysis["experimental_setups"]):
            LOGGER.info("{0:*^78s}".format(cntrl_type))
            namespace["prepared"][version][cntrl_type] = dict()
            for (exp_name, exp_setup) in izip(experiments, setups):
                LOGGER.info("{0:*^78s}".format(exp_name))
                df = organism.activity[exp_name]
                setup_func = glob_vars[exp_setup]
                namespace["prepared"][version][cntrl_type][exp_name] =\
                        setup_func(cntrl_type, df, id2gene)
                if any(method.startswith("delayed") for method in chain(analysis["control"], *analysis["ctc"])):
                    namespace["prepared"][version][cntrl_type][exp_name]["delayed"] = dict()
                    for delta in analysis["delays"]:
                        delayed_continuous(namespace["prepared"][version][cntrl_type][exp_name],
                                delta)
                if any(ms_name.endswith("comparison") for ms_name in chain(*analysis["measures"])):
                    rate_continuous(namespace["prepared"][version][cntrl_type][exp_name])
    LOGGER.debug("\n".join(print_dict(namespace)))
    db = shelve.open(config["shelve"], protocol=pickle.HIGHEST_PROTOCOL)
    for (key, value) in namespace.iteritems():
        db[key] = value
    db.close()
    # general parallel setup using IPython.parallel
    LOGGER.info("Remote imports")
    d_view = remote_client.direct_view()
    d_view.execute("import numpy as np; "\
            "import shelve; import pickle;"\
            "import pyorganism as pyorg; import pyorganism.regulation as pyreg;"\
            "import logging; from IPython.config import Application;"\
            "LOGGER = Application.instance().log;"\
            "LOGGER.setLevel(logging.{level});".format(level=args.log_level),
            block=True)
    LOGGER.info("Transfer data")
#    d_view.push(namespace, block=True)
    d_view.execute("db = shelve.open('{shelve}', protocol=pickle.HIGHEST_PROTOCOL);"\
                    "globals().update(db);db.close()".format(shelve=config["shelve"]), block=True)
    LOGGER.info("Generate job descriptions")
    jobs = job_gen(organism, config, namespace)
    l_view = remote_client.load_balanced_view()
    bar = ProgressBar(maxval=len(jobs), widgets=[Timer(), " ", Percentage(),
            " ", Bar(), " ", ETA()]).start()
    result_mngr = ResultManager(config["output"], table_key)
    results_it = l_view.map(worker, jobs, ordered=False, block=False)
    for (spec, res_cntrl, res_ctc, samples) in results_it:
        LOGGER.debug(res_cntrl)
        LOGGER.debug(res_ctc)
        result(result_mngr, spec, res_cntrl, res_ctc, samples)
        bar += 1
    result_mngr.finalize()
    bar.finish()
    LOGGER.info("parallel speed-up was %.3g",
            results_it.serial_time / results_it.wall_time)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=None)
    parser.add_argument("-v", "--version", action="version", version="0.1")
    parser.add_argument("--config", dest="config", default="config.json",
            help="Configuration file to use (default: %(default)s)")
    parser.add_argument("--profile", dest="profile", default="default",
            help="IPython profile to connect to cluster (default: %(default)s)")
    parser.add_argument("--cluster-id", dest="cluster_id", default=None,
            help="IPython cluster-id to connect to (default: %(default)s)")
    parser.add_argument("--log-level", dest="log_level", default="INFO",
            help="Log level, i.e., DEBUG, INFO, WARN, ERROR, CRITICAL (default: %(default)s)")
    parser.add_argument("--encoding", dest="encoding", default="utf-8",
            help="File encoding to assume (default: %(default)s)")
    args = parser.parse_args()
    dictConfig({"version": 1, "incremental": True, "root": {"level": args.log_level}})
    remote_client = Client(profile=args.profile, cluster_id=args.cluster_id)
    pid_map = remote_client[:].apply_async(os.getpid).get_dict()
    LOGGER.debug(str(pid_map))
    try:
        sys.exit(main(remote_client, args))
    except: # we want to catch everything
        (err, msg, trace) = sys.exc_info()
        if not isinstance(err, (SystemExit, RemoteError)):
            for (engine_id, pid) in pid_map.iteritems():
                LOGGER.debug("interrupting engine %d", engine_id)
                try:
                    os.kill(pid, SIGINT)
                except OSError:
                    continue
            raise err, msg, trace
        if isinstance(err, RemoteError):
            raise err, msg, trace
    finally:
        logging.shutdown()

