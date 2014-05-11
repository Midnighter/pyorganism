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

from itertools import izip

from IPython.parallel import (Client, interactive)
from progressbar import ProgressBar


logging.basicConfig()
LOGGER = logging.getLogger()
#LOGGER.setLevel(logging.INFO)
LOGGER.setLevel(logging.DEBUG)
#LOGGER.addHandler(logging.StreamHandler())


def check_path(path):
    if not os.path.exists(path):
        raise IOError("file does not exist '{path}'".format(path=path))

def simple_continuous(df, feature2gene):
    df[df <= 0.0] = numpy.nan
    df = df.apply(pyorg.norm_zero2unity, axis=1, raw=True)
    results = dict()
    results["time"] = list()
    results["active"] = list()
    results["levels"] = list()
    for col in df.columns:
        eligible = df[col][numpy.isfinite(df[col])]
        mask = [bool(feature2gene[name]) for name in eligible.index]
        active = [feature2gene[name] for name in eligible.index[mask]]
        levels = eligible[mask]
        LOGGER.info("        %s min: %d active genes", col, len(active))
        results["time"].append(int(float(col)))
        results["active"].append(active)
        results["levels"].append(levels)
    return results

def shuffle(df, n_times=1, axis=0):
    df = df.copy()
    axis ^= 1 # make axes between numpy and pandas behave as expected
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
    for (basename, name, reader, extra) in izip(config["experimental_paths"],
            config["experimental_names"], config["experimental_reader"],
            config["experimental_args"]):
        reader_func = getattr(pymicro, reader)
        path = os.path.join(config["experimental_base"], basename)
        check_path(path)
        LOGGER.info("  %s: '%s'", name, path)
        organism.activity[name] = reader_func(path, mutants=extra)

def continuous_jobs(organism, config):
    LOGGER.info("Generating continuous job specifications:")
    jobs = list()
    for version in config["versions"]:
        for (cntrl_name, gene_file, map_file, net_file, experiments, setups,
                method, measures, random_num, robustness_num,
                rob_extra) in izip(config["control_types"],
                config["gene_paths"], config["mapping_paths"], config["network_paths"],
                config["experimental_sets"], config["experimental_setup"], config["method"],
                config["measures"], config["random_num"],
                config["robustness_num"], config["robustness_args"]):
            gene_path = os.path.join(config["regulondb_base"], version, gene_file)
            check_path(gene_path)
            map_path = os.path.join(config["regulondb_base"], version, map_file)
            check_path(map_path)
            net_path = os.path.join(config["regulondb_base"], version, net_file)
            check_path(net_path)
            for ms_name in measures:
                for (exp_name, exp_setup) in izip(experiments, setups):
                    times = organism.activity[exp_name].columns
                    for time_point in times:
                        spec = dict()
                        spec["version"] = version
                        spec["control"] = cntrl_name
                        spec["time"] = time_point
                        spec["genes"] = gene_path
                        spec["mapping"] = map_path
                        spec["network"] = net_path
                        spec["experiment"] = exp_name
                        spec["setup"] = exp_setup
                        spec["method"] = method
                        spec["measure"] = ms_name
                        spec["random_num"] = random_num
                        spec["robustness_num"] = robustness_num
                        spec["robustness_args"] = rob_extra
                        jobs.append(spec)
    LOGGER.info("  %d jobs", len(jobs))
    return jobs

@interactive
def worker(spec):
    organism = globals()["organism"]
    setup_func = globals()[spec["setup"]]
    method = getattr(pyreg, spec["method"])
    measure = getattr(pyreg, spec["measure"])
    genes = pyorg.read_pickle(spec["genes"]) # load genes into memory
    id2gene = pyorg.read_pickle(spec["mapping"])
    net = pyorg.read_pickle(spec["network"])
    df = organism.activity[spec["experiment"]]
    prepared = setup_func(df, id2gene)
    index = prepared["time"].index(spec["time"])
    results = method(net, prepared["active"][index], prepared["levels"][index],
            random_num=spec["random_num"], measure=measure)
    pyreg.clear_memory()
    return results

def main(args):
    config = json.load(codecs.open(args.config, encoding=args.encoding, mode="rb"))
    organism = pyorg.Organism(name=config["organism"])
    if config["continuous"]:
        load_continuous(organism, config)
    else:
        load_discrete(organism, config)
    # general parallel setup using IPython.parallel
    rc = Client(profile=args.profile, cluster_id=args.cluster_id)
    dv = rc.direct_view()
    LOGGER.info("Remote imports")
    dv.execute("import numpy; "\
            "import pyorganism as pyorg; import pyorganism.regulation as pyreg;"\
            "import logging; from IPython.config import Application;"\
            "LOGGER = Application.instance().log;"\
            "LOGGER.setLevel(logging.DEBUG);", block=True) # make level variable
    LOGGER.info("Transfer data")
    if config["continuous"]:
        dv.push({"organism": organism, "simple_continuous": simple_continuous,
                "randomised_continuous": randomised_continuous,
                "fully_randomised_continuous": fully_randomised_continuous},
                block=True)
    else:
        dv.push({"organism": organism}, block=True)
    if config["continuous"]:
        jobs = continuous_jobs(organism, config)
    else:
        jobs = discrete_jobs(organism, config)
    lv = rc.load_balanced_view()
    bar = ProgressBar(maxval=len(jobs), widgets=None).start()
    results_it = lv.map(worker, jobs, ordered=False, block=False)
    for result in results_it:
        bar.update(bar.currval + 1)
    bar.finish()
    LOGGER.info("parallel speed-up was %.3g",
            results_it.serial_time / results_it.wall_time)
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=None)
    parser.add_argument("-v", "--version", action="version", version="0.1")
    parser.add_argument("--config", dest="config", default="config.json",
            help="Configuration file to use (default: %(default)s)")
    parser.add_argument("--profile", dest="profile", default="default",
            help="IPython profile to connect to cluster (default: %(default)s)")
    parser.add_argument("--cluster-id", dest="cluster_id", default=None,
            help="IPython cluster-id to connect to (default: %(default)s)")
    parser.add_argument("--level", dest="level", default="info",
            help="Log level, i.e., debug, info, warn, error, critical (default: %(default)s)")
    parser.add_argument("--encoding", dest="encoding", default="utf-8",
            help="File encoding to assume (default: %(default)s)")
    args = parser.parse_args()
    sys.exit(main(args))

