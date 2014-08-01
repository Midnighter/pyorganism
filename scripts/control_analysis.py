# -*- coding: utf-8 -*-


import sys
import os
import logging
import argparse
import json
import codecs

import numpy

import pyorganism as pyorg
import pyorganism.regulation as pyreg
import pyorganism.io.microarray as pymicro

from itertools import (izip,)
from signal import SIGINT
from exceptions import SystemExit
from logging.config import dictConfig

from IPython.parallel import (Client, interactive)
from progressbar import (ProgressBar, Timer, Bar, Percentage, ETA)

from pyorganism.io.hdf5 import ResultManager


LOGGER = logging.getLogger()
LOGGER.addHandler(logging.StreamHandler())


def check_path(path):
    if not os.path.exists(path):
        raise IOError("file does not exist '{path}'".format(path=path))


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

def simple_discrete(df, name2gene):
    eligible = df["name"].notnull()
    active = [name2gene[name] for name in df["name"][eligible]\
            if name2gene[name] is not None]
    LOGGER.info("  mapped %d/%d active genes", len(active), len(df))
    return active

def ratio_discrete(df, name2gene, direction="up"):
    eligible = df["name"].notnull()
    if direction == "up":
        mask = (df["ratio (A/B)"] > 1.0) & eligible
    elif direction == "down":
        mask = (df["ratio (A/B)"] < 1.0) & eligible
    else:
        mask = eligible
    active = [name2gene[name] for name in df["name"][mask]\
            if name2gene[name] is not None]
    LOGGER.info("  mapped %d/%d {dir}-regulated genes", len(active), len(df))
    return active

def discrete_jobs(organism, config):
    LOGGER.info("Generating discrete job specifications:")
    jobs = list()
    analysis = config["analysis"]
    for version in config["versions"]:
        for (cntrl_name, experiments, setups,
                control, ctc, measures, random_num, robustness_num,
                rob_extra) in izip(analysis["control_types"],
                analysis["experimental_sets"], analysis["experimental_setups"],
                analysis["control"], analysis["ctc"],
                analysis["measures"], analysis["random_num"],
                analysis["robustness_num"], analysis["robustness_args"]):
            if cntrl_name == "digital":
                index = config["network"]["names"].index("TRN")
            elif cntrl_name == "analog":
                index = config["network"]["names"].index("GPN")
            projections = config["network"]["projections"][index]
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
                                spec["basis"] = basis
                                spec["direction"] = direction
                                spec["control"] = control
                                spec["ctc"] = ctc
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
                            spec["basis"] = basis
                            spec["control"] = control
                            spec["ctc"] = ctc
                            spec["measure"] = ms_name
                            spec["random_num"] = random_num
                            spec["robustness_num"] = robustness_num
                            spec["robustness_args"] = rob_extra
                            jobs.append(spec)
    LOGGER.info("  %d jobs", len(jobs))
    return jobs

@interactive
def discrete_worker(spec):
    version = spec["version"]
    cntrl_type = spec["control_type"]
    global_vars = globals()
    control = getattr(pyreg, spec["control"])
    ctc = getattr(pyreg, spec["ctc"])
    measure = getattr(pyreg, spec["measure"])
    net = global_vars["networks"][version][cntrl_type][spec["projection"]]
    if cntrl_type == "analog":
        active = global_vars["prepared"][version][cntrl_type][spec["experiment"]][spec["direction"]]
    else:
        active = global_vars["prepared"][version][cntrl_type][spec["experiment"]]
    LOGGER.debug(len(active))
    res_cntrl = control(net, active, measure=measure)
    res_ctc = ctc(net, active, random_num=spec["random_num"], measure=measure)
    pyreg.clear_memory()
    return (spec, res_cntrl, res_ctc)

def discrete_result(manager, spec, res_cntrl, res_ctc):
    manager.append(spec["version"], spec["control_type"], spec["continuous"],
            spec["experiment"], res_cntrl, spec["control"], res_ctc,
            spec["ctc"], spec["measure"], direction=spec.get("direction", None))


##############################################################################
# Continuous
##############################################################################


def simple_continuous(df, feature2gene):
    df[df <= 0.0] = numpy.nan
    df = df.apply(pyorg.norm_zero2unity, axis=1, raw=True)
    results = dict()
    results["time"] = list()
    results["active"] = list()
    results["levels"] = list()
    for col in df.columns:
        eligible = df[col][numpy.isfinite(df[col])]
        mask = [feature2gene[name] is not None for name in eligible.index]
        active = [feature2gene[name] for name in eligible[mask].index]
        levels = eligible[mask]
        LOGGER.info("        %s min: %d active genes", col, len(active))
        results["time"].append(col)
        results["active"].append(active)
        results["levels"].append(levels)
    return results

def rate_continuous(df, feature2gene):
    df[df <= 0.0] = numpy.nan
    df = df.apply(pyorg.norm_zero2unity, axis=1, raw=True)
    results = dict()
    results["time"] = list()
    results["active"] = list()
    results["levels"] = list()
    num_cols = len(df.columns)
    for i in range(num_cols - 1):
        col_a = df.icol(i)
        col_b = df.icol(i + 1)
        # TODO: fix selection of names and values
#        eligible = df.index[numpy.isfinite(col_a) & numpy.isfinite(col_b)]
#        mask = [feature2gene[name] is not None for name in eligible.index]
#        active = [feature2gene[name] for name in eligible[mask].index]
#        levels = col_b[eligible[mask]] - col_a.loc[eligible[mask]]
        LOGGER.info("        %s - %s min: %d active genes", col_a.name,
                col_b.name, len(active))
        results["time"].append(col_b.name)
        results["active"].append(active)
        results["levels"].append(levels)
    return results

def shuffle(df, n_times=1, axis=0):
    df = df.copy()
    axis ^= 1 # make axes between numpy and pandas behave as expected (only 2D!)
    for _ in range(n_times):
        for view in numpy.rollaxis(df.values, axis):
            numpy.random.shuffle(view)
    return df

def randomised_continuous(df, feature2gene):
    df = shuffle(df)
    return simple_continuous(df, feature2gene)

def fully_randomised_continuous(df, feature2gene):
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

def continuous_jobs(organism, config):
    LOGGER.info("Generating continuous job specifications:")
    jobs = list()
    analysis = config["analysis"]
    for version in config["versions"]:
        for (cntrl_name, experiments, control, ctc, measures, random_num,
                robustness_num, rob_extra) in izip(analysis["control_types"],
                analysis["experimental_sets"], analysis["control"],
                analysis["ctc"], analysis["measures"], analysis["random_num"],
                analysis["robustness_num"], analysis["robustness_args"]):
            if cntrl_name == "digital":
                index = config["network"]["names"].index("TRN")
            elif cntrl_name == "analog":
                index = config["network"]["names"].index("GPN")
            projections = config["network"]["projections"][index]
            for basis in projections:
                for ms_name in measures:
                    for exp_name in experiments:
    # TODO: need to change this to use the actual times (for rate and delayed)
                        times = organism.activity[exp_name].columns
                        for time_point in times:
                            spec = dict()
                            spec["version"] = version
                            spec["continuous"] = config["continuous"]
                            spec["control_type"] = cntrl_name
                            spec["time"] = time_point
                            spec["experiment"] = exp_name
                            spec["basis"] = basis
                            spec["control"] = control
                            spec["ctc"] = ctc
                            spec["measure"] = ms_name
                            spec["random_num"] = random_num
                            spec["robustness_num"] = robustness_num
                            spec["robustness_args"] = rob_extra
                            jobs.append(spec)
    LOGGER.info("  %d jobs", len(jobs))
    return jobs

@interactive
def continuous_worker(spec):
    version = spec["version"]
    cntrl_type = spec["control_type"]
    global_vars = globals()
    control = getattr(pyreg, spec["control"])
    ctc = getattr(pyreg, spec["ctc"])
    measure = getattr(pyreg, spec["measure"])
    # TODO: projection
    net = global_vars["networks"][version][cntrl_type]
    prepared = global_vars["prepared"][version][cntrl_type][spec["experiment"]]
    index = prepared["time"].index(spec["time"])
    active = prepared["active"][index]
    LOGGER.debug(len(active))
    levels = prepared["levels"][index]
    res_cntrl = control(net, active, levels, measure=measure)
    res_ctc = ctc(net, active, levels, random_num=spec["random_num"], measure=measure)
    return (spec, res_cntrl, res_ctc)

def continuous_result(manager, spec, res_cntrl, res_ctc):
    manager.append(spec["version"], spec["control_type"], spec["continuous"],
            spec["experiment"], res_cntrl, spec["control"], res_ctc,
            spec["ctc"], spec["measure"], time=int(float(spec["time"])))


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
    namespace["id2gene"] = dict()
    namespace["networks"] = dict()
    namespace["prepared"] = dict()
    for version in config["versions"]:
        namespace["genes"][version] = pyorg.read_pickle(os.path.join(
                data["base"], version, data["gene_path"]))
        id2gene = pyorg.read_pickle(os.path.join(data["base"], version,
                data["mapping_path"]))
        namespace["networks"][version] = dict()
        for (net_name, net_file, projections) in izip(network["names"],
                network["paths"], network["projections"]):
            net = pyorg.read_pickle(os.path.join(data["base"], version, net_file))
            namespace["networks"][version][net_name] = dict()
            for basis in projections:
                if basis == "gene":
                    namespace["networks"][version][net_name][basis] = net
                elif basis == "transcription_unit":
                    namespace["networks"][version][net_name][basis] =\
                            pyreg.to_transcription_unit_based(net)
                elif basis == "operon":
                    namespace["networks"][version][net_name][basis] =\
                            pyreg.to_operon_based(net)
        namespace["prepared"][version] = dict()
        for (cntrl_type, experiments, setups) in izip(analysis["control_types"],
                 analysis["experimental_sets"], analysis["experimental_setup"]):
            namespace["prepared"][version][cntrl_type] = dict()
            for (exp_name, exp_setup) in izip(experiments, setups):
                df = organism.activity[exp_name]
                setup_func = glob_vars[exp_setup]
                if exp_setup == "ratio_discrete":
                    namespace["prepared"][version][cntrl_type][exp_name] = dict()
                    namespace["prepared"][version][cntrl_type][exp_name]["up"] =\
                            setup_func(df, id2gene, "up")
                    namespace["prepared"][version][cntrl_type][exp_name]["down"] =\
                            setup_func(df, id2gene, "down")
                else:
                    namespace["prepared"][version][cntrl_type][exp_name] =\
                            setup_func(df, id2gene)
    # general parallel setup using IPython.parallel
    LOGGER.info("Remote imports")
    d_view = remote_client.direct_view()
    d_view.execute("import numpy; "\
            "import pyorganism as pyorg; import pyorganism.regulation as pyreg;"\
            "import logging; from IPython.config import Application;"\
            "LOGGER = Application.instance().log;"\
            "LOGGER.setLevel(logging.{level});".format(level=args.log_level),
            block=True)
    LOGGER.info("Transfer data")
    d_view.push(namespace, block=True)
    result_mngr = ResultManager(config["output"], table_key)
    jobs = job_gen(organism, config)
    l_view = remote_client.load_balanced_view()
    bar = ProgressBar(maxval=len(jobs), widgets=[Timer(), " ", Percentage(),
            " ", Bar(), " ", ETA()]).start()
    results_it = l_view.map(worker, jobs, ordered=False, block=False)
    for (spec, res_cntrl, res_ctc) in results_it:
        LOGGER.debug(res_cntrl)
        LOGGER.debug(res_ctc)
        result(result_mngr, spec, res_cntrl, res_ctc)
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
        if not isinstance(err, SystemExit):
            for (engine_id, pid) in pid_map.iteritems():
                LOGGER.debug("interrupting engine %d", engine_id)
                try:
                    os.kill(pid, SIGINT)
                except OSError:
                    continue
            raise err, msg, trace
    finally:
        logging.shutdown()

