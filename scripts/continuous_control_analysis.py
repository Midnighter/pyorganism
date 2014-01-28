# -*- coding: utf-8 -*-


"""
================
Control Analysis
================

:Authors:
    Moritz Emanuel Beber
:Date:
    2012-07-30
:Copyright:
    Copyright(c) 2012 Jacobs University of Bremen. All rights reserved.
:File:
    control_analysis.py
"""


import sys
import os
import logging

import numpy

import pyorganism
import pyorganism.regulation as pyreg
import pyorganism.regulation.parallel as pp

from itertools import izip

from IPython.parallel import Client

from pyorganism.io.hdf5 import ResultManager
from pyorganism.io.microarray import read_interpolated
from pyorganism.statistics import norm_zero_unity


LOGGER = logging.getLogger()
LOGGER.addHandler(logging.StreamHandler())
LOGGER.setLevel(logging.INFO)


################################################################################
# functions known at configuration time

def simple_continuous(df, feature2gene):
    df[df <= 0.0] = numpy.nan
    df = df.apply(norm_zero_unity, axis=1, raw=True)
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

def shuffle(df, axis=0):
    df = df.copy()
    df.apply(numpy.random.shuffle, axis=axis, raw=False)
    return df

def randomised_continuous(df, feature2gene):
    df = shuffle(df)
    return simple_continuous(df, feature2gene)

def fully_randomised_continuous(df, feature2gene):
    df = shuffle(df)
    df = shuffle(df, axis=1)
    return simple_continuous(df, feature2gene)

################################################################################
# this should be moved to JSON
base_dir = "Expression/intra_strain"
CONFIG = dict(
        continuous=True,
        data_paths=[
            (os.path.join(base_dir, "time_course_wt_interpol.txt"), "wt"),
            (os.path.join(base_dir, "time_course_fis_interpol.txt"), "fis"),
            (os.path.join(base_dir, "time_course_hns_interpol.txt"), "hns")
        ],
        data_load=[
            read_interpolated,
            read_interpolated,
            read_interpolated
        ],
        data_args=[
            [[]],
            [["fis"]],
            [["hns"]]
        ],
# the type of analysis to perform on the data
        analyses=["digital", "analog"],
# choice of data sets per analysis
        data_sets=[
                ['wt', 'fis', 'hns'],
                ['wt', 'fis', 'hns']
        ],
        significance=[
                simple_continuous,
                simple_continuous
        ],
        evaluater=[
            [pyreg.continuous_difference_coherence,
                pyreg.continuous_abs_coherence,
                pyreg.continuous_functional_coherence,
                pyreg.continuous_functional_comparison],
            [pyreg.continuous_difference_coherence,
                pyreg.continuous_abs_coherence]
        ],
# other parameters, per analysis
        random_num=[
                1E5,
                1E5
        ],
        jackknife_num=[
                0,
                0
        ],
        jackknife_fraction=[
                0.1,
                0.1
        ]
)
################################################################################

def prepare_continuous_digital(db_path, organism, d_view):
    LOGGER.info("Loading TRN")
    organism.trn = pyorganism.read_pickle(os.path.join(db_path, "trn.pkl"))
    num = sum(1 for node in organism.trn.nodes_iter() if isinstance(node,
            pyreg.TranscriptionFactor))
    LOGGER.info("%d transcription factors (%3.2f%%)", num, 100.0 * num /
            len(organism.trn))
    num = len(organism.trn) - num
    LOGGER.info("%d regulated genes (%3.2f%%)", num, 100.0 * num /
            len(organism.trn))
    LOGGER.info("%d regulatory links", organism.trn.size())
    pp.prepare_continuous_digital(d_view)

def prepare_continuous_analog(db_path, organism, d_view):
    LOGGER.info("Loading GPN")
    organism.gpn = pyorganism.read_pickle(os.path.join(db_path, "gpn_5000.pkl"))
    LOGGER.info("%d genes", len(organism.gpn))
    LOGGER.info("%d links", organism.gpn.size())
    pp.prepare_continuous_analog(d_view)

CONTROL_FUNCS = dict(
        digital=pyreg.continuous_digital_control,
        analog=pyreg.continuous_analog_control
)

CTC_FUNCS = dict(
        digital=pp.continuous_digital_ctc,
        analog=pp.continuous_analog_ctc
)

def anonymous_continuous(name, version, data_name, network, active, levels,
        random_num, jackknife_fraction, jackknife_num, evaluate, time,
        results, dv, lv=None, note=""):
    LOGGER.info("Computing %s control", name)
    strength = CONTROL_FUNCS[name](network, active, levels)
    LOGGER.info("Computing %s ctc", name)
    ctc = CTC_FUNCS[name](dv, network, active, levels, random_num, lb_view=lv)
    if jackknife_num > 0:
        LOGGER.info("Computing %s ctc robustness", name)
        LOGGER.setLevel(logging.WARN)
        jack_num = int(numpy.floor(len(active) * jackknife_fraction))
        robust = list()
        for i in range(jackknife_num):
            jacked_ind = range(len(active))
            jacked_ind = [jacked_ind.pop(numpy.random.randint(len(jacked_ind))) for j in range(jack_num)]
            jacked = [active[j] for j in jacked_ind]
            jacked_levels = [levels[j] for j in jacked_ind]
            z_score = CTC_FUNCS[name](dv, network, jacked, jacked_levels, random_num, lb_view=lv)
            robust.append(z_score)
        LOGGER.setLevel(logging.INFO)
        robust = numpy.array(robust, dtype=float)
        results.append(version, name, False, data_name, strength, ctc,
                robustness=robust, time=time, note=note)
    else:
        results.append(version, name, False, data_name, strength, ctc, time=time, note=note)

def digital_continuous(name, version, data_name, organism, random_num, jackknife_fraction,
        jackknife_num, evaluaters, results, dv, lv=None, note=""):
    LOGGER.info("*" * 79)
    LOGGER.info(data_name)
    times = organism.significant[name][data_name]["time"]
    active_genes = organism.significant[name][data_name]["active"]
    levels = organism.significant[name][data_name]["levels"]
    for eval_func in evaluaters:
        LOGGER.info("*" * 79)
        LOGGER.info("Evaluater: %s", eval_func.__name__)
        for i in range(len(active_genes)):
            LOGGER.info("*" * 79)
            LOGGER.info("Version: %s Time: %d", version, times[i])
            LOGGER.info("*" * 79)
            anonymous_continuous(name, version, data_name, organism.trn,
                active_genes[i], levels[i], random_num, jackknife_fraction,
                jackknife_num, eval_func, times[i], results, dv, lv,
                note=eval_func.__name__)

def analog_continuous(name, version, data_name, organism, random_num, jackknife_fraction,
        jackknife_num, evaluaters, results, dv, lv=None, note=""):
    LOGGER.info("*" * 79)
    LOGGER.info(data_name)
    times = organism.significant[name][data_name]["time"]
    active_genes = organism.significant[name][data_name]["active"]
    levels = organism.significant[name][data_name]["levels"]
    for eval_func in evaluaters:
        LOGGER.info("*" * 79)
        LOGGER.info("Evaluater: %s", eval_func.__name__)
        for i in range(len(active_genes)):
            LOGGER.info("*" * 79)
            LOGGER.info("Version: %s Time: %d", version, times[i])
            LOGGER.info("*" * 79)
            anonymous_continuous(name, version, data_name, organism.gpn,
                active_genes[i], levels[i], random_num, jackknife_fraction,
                jackknife_num, eval_func, times[i], results, dv, lv,
                note=eval_func.__name__)

ANAL_PREP = dict(
        digital=prepare_continuous_digital,
        analog=prepare_continuous_analog,
#        metabolic=prepare_metabolic
)

ANAL_RESOLVE = dict(
        digital=digital_continuous,
        analog=analog_continuous,
#        metabolic=prepare_metabolic
)

def load_continuous(db_path, organism, experiment_paths, loading_funcs,
        data_args):
    LOGGER.info("Loading continuous expression data:")
    for ((path, name), load_data, extra_args) in izip(experiment_paths,
            loading_funcs, data_args):
        LOGGER.info("    %s: '%s'", name, path)
        organism.activity[name] = load_data(path, *extra_args)

def prepare_continuous(db_path, organism, analyses, analysis_sets, significance_functions):
    LOGGER.info("Loading genes")
    organism.genes = pyorganism.read_pickle(os.path.join(db_path, "genes.pkl"))
    LOGGER.info("Found %d", len(organism.genes))
    LOGGER.info("Loading feature to gene map")
    feature2gene = pyorganism.read_pickle(os.path.join(db_path, "feature2gene.pkl"))
    num = sum(1 for gene in feature2gene.itervalues() if gene)
    LOGGER.info("Map contains %d features and %d genes (%3.2f%%)", len(feature2gene),
            num, 100.0 * num / len(feature2gene))
    for (analysis, data_sets, sig_func) in izip(analyses, analysis_sets, significance_functions):
        LOGGER.info("Analysis '%s':", analysis)
        organism.significant[analysis] = dict()
        for name in data_sets:
            LOGGER.info("    Experiment '%s':", name)
            df = organism.activity[name]
            organism.significant[analysis][name] = sig_func(df, feature2gene)

def main(db_path, output):
    (base_dir, version) = os.path.split(db_path)
    if not version:
        (base_dir, version) = os.path.split(os.path.dirname(db_path))
    ecoli = pyorganism.Organism(name="E. coli K12")
    if CONFIG["continuous"]:
        load_continuous(db_path, ecoli, CONFIG["data_paths"],
                CONFIG["data_load"], CONFIG["data_args"])
        prepare_continuous(db_path, ecoli, CONFIG["analyses"], CONFIG["data_sets"],
                CONFIG["significance"])
        results = ResultManager(output, "/Continuous")
    else:
        pass
    # general parallel setup using IPython.parallel
    rc = Client()
    dv = rc.direct_view()
#    lv = rc.load_balanced_view()
    lv = None
    LOGGER.info("Remote imports")
    dv.execute("import pyorganism;import pyorganism.regulation as pyreg", block=True)
    for (analysis, data, random_num, jackknife_fraction, jackknife_num, evaluater) in izip(
            CONFIG["analyses"], CONFIG["data_sets"], CONFIG["random_num"],
            CONFIG["jackknife_fraction"], CONFIG["jackknife_num"], CONFIG["evaluater"]):
        LOGGER.info("*" * 79)
        LOGGER.info(analysis.upper())
        LOGGER.info("*" * 79)
        analysis_function = ANAL_RESOLVE.get(analysis, StandardError("unknown name"))
        ANAL_PREP.get(analysis, StandardError("unknown name"))(db_path, ecoli, dv)
        for name in data:
            analysis_function(analysis, version, name, ecoli, random_num,
                    jackknife_fraction, jackknife_num, evaluater, results, dv, lv)
    results.finalize()


if __name__ == "__main__":
    if len(sys.argv) != 3:
        LOGGER.critical("%s <RegulonDB objects: path> <results output: path>", sys.argv[0])
        sys.exit(2)
    else:
        main(sys.argv[1], sys.argv[2])


