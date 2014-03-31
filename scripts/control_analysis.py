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

from copy import copy
from itertools import izip

from IPython.parallel import Client

from pyorganism.io.hdf5 import ResultManager
from pyorganism.io.microarray import read_microarray


LOGGER = logging.getLogger()
LOGGER.addHandler(logging.StreamHandler())
LOGGER.setLevel(logging.INFO)


################################################################################
# functions known at configuration time

def simple_discrete(df, name2gene, blattner2gene):
    active = list()
    for (name, blattner) in df[["name", "blattner"]].itertuples(index=False):
        name_gene = name2gene.get(name)
        if name_gene:
            active.append(name_gene)
            continue
        blattner_gene = blattner2gene.get(blattner)
        if blattner_gene:
            active.append(blattner_gene)
            continue
    LOGGER.info("        %d DEGs", len(active))
    return active

def high_low_ratio(df, name2gene, blattner2gene):
    pos = list()
    neg = list()
    for (name, blattner, ratio) in df[["name", "blattner", "ratio (A/B)"]].itertuples(index=False):
        name_gene = name2gene.get(name)
        if name_gene:
            if ratio > 1.0:
                pos.append(name_gene)
            if ratio < 1.0:
                neg.append(name_gene)
            continue
        blattner_gene = blattner2gene.get(blattner)
        if blattner_gene:
            if ratio > 1.0:
                pos.append(blattner_gene)
            if ratio < 1.0:
                neg.append(blattner_gene)
            continue
    LOGGER.info("        %d over-expressed genes", len(pos))
    LOGGER.info("        %d under-expressed genes", len(neg))
    return (pos, neg)

################################################################################
# this should be moved to JSON
base_data = "Expression/LZ41-LZ54_single_knockouts"
CONFIG = dict(
        continuous=False,
        data_paths=[
            (os.path.join(base_data, "LZ41-LZ54.tsv"), "wt"),
            (os.path.join(base_data, "LZ41_d_fis-LZ54_d_fis.tsv"), "fis"),
            (os.path.join(base_data, "LZ41_d_hns-LZ54_d_hns.tsv"), "hns"),
            (os.path.join(base_data, "LZ41-LZ41_d_fis.tsv"), "wt-fis-low"),
            (os.path.join(base_data, "LZ54-LZ54_d_fis.tsv"), "wt-fis-high"),
            (os.path.join(base_data, "LZ41-LZ41_d_hns.tsv"), "wt-hns-low"),
            (os.path.join(base_data, "LZ54-LZ54_d_hns.tsv"), "wt-hns-high")
        ],
        data_load=[
            read_microarray,
            read_microarray,
            read_microarray,
            read_microarray,
            read_microarray,
            read_microarray,
            read_microarray
        ],
# the type of analysis to perform on the data
        analyses=["digital", "analog"],
# choice of data sets per analysis
        data_sets=[
                ['wt', 'fis', 'hns', 'wt-fis-low', 'wt-fis-high', 'wt-hns-low',
                        'wt-hns-high'],
                ['wt-fis-low', 'wt-fis-high', 'wt-hns-low', 'wt-hns-high']
        ],
        significance=[
                simple_discrete,
                high_low_ratio
        ],
# other parameters, per analysis
        random_num=[
                1E6,
                1E6
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

def prepare_digital(db_path, organism, d_view):
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
    pp.prepare_digital(d_view, organism.trn)

def prepare_analog(db_path, organism, d_view):
    LOGGER.info("Loading GPN")
    organism.gpn = pyorganism.read_pickle(os.path.join(db_path, "gpn_5000.pkl"))
    LOGGER.info("%d genes", len(organism.gpn))
    LOGGER.info("%d links", organism.gpn.size())
    pp.prepare_analog(d_view, organism.gpn)

CONTROL_FUNCS = dict(
        digital=pyreg.digital_control,
        analog=pyreg.analog_control
)

CTC_FUNCS = dict(
        digital=pp.digital_ctc,
        analog=pp.analog_ctc
)

def anonymous_discrete(name, version, data_name, network, active, random_num, jackknife_fraction,
        jackknife_num, results, dv, lv=None, note=""):
    LOGGER.info("*" * 79)
    LOGGER.info(data_name)
    LOGGER.info("*" * 79)
    LOGGER.info("Computing %s control", name)
    strength = CONTROL_FUNCS[name](network, active)
    LOGGER.info("Computing %s ctc", name)
    ctc = CTC_FUNCS[name](dv, network, active, random_num, lb_view=lv)
    if jackknife_num > 0:
        LOGGER.info("Computing %s ctc robustness", name)
        LOGGER.setLevel(logging.WARN)
        jack_num = int(numpy.floor(len(active) * jackknife_fraction))
        robust = list()
        for i in range(jackknife_num):
            jacked = copy(active)
            jacked = [jacked.pop(numpy.random.randint(len(jacked))) for j in
                range(jack_num)]
            z_score = CTC_FUNCS[name](dv, network, active, random_num, lb_view=lv)
            robust.append(z_score)
        LOGGER.setLevel(logging.INFO)
        robust = numpy.array(robust, dtype=float)
        results.append(version, name, False, data_name, strength, ctc,
                robustness=robust, note=note)
    else:
        results.append(version, name, False, data_name, strength, ctc, note=note)

def digital_discrete(name, version, data_name, organism, random_num, jackknife_fraction,
        jackknife_num, results, dv, lv=None, note=""):
    anonymous_discrete(name, version, data_name, organism.trn,
            organism.significant[name][data_name], random_num, jackknife_fraction,
            jackknife_num, results, dv, lv, note="")

def analog_discrete(name, version, data_name, organism, random_num, jackknife_fraction,
        jackknife_num, results, dv, lv=None, note=""):
    (pos, neg) = organism.significant[name][data_name]
    anonymous_discrete(name, version, data_name, organism.gpn, pos, random_num,
            jackknife_fraction, jackknife_num, results, dv, lv, note="over")
    anonymous_discrete(name, version, data_name, organism.gpn, neg, random_num,
            jackknife_fraction, jackknife_num, results, dv, lv, note="under")

ANAL_PREP = dict(
        digital=prepare_digital,
        analog=prepare_analog,
#        metabolic=prepare_metabolic
)

ANAL_RESOLVE = dict(
        digital=digital_discrete,
        analog=analog_discrete,
#        metabolic=prepare_metabolic
)

def load_discrete(db_path, organism, experiment_paths, loading_funcs):
    LOGGER.info("Loading discrete expression data:")
    for ((path, name), load_data) in izip(experiment_paths, loading_funcs):
        LOGGER.info("    %s: '%s'", name, path)
        organism.activity[name] = load_data(path)

def prepare_discrete(db_path, organism, analyses, analysis_sets, significance_functions):
    LOGGER.info("Loading genes")
    organism.genes = pyorganism.read_pickle(os.path.join(db_path, "genes.pkl"))
    LOGGER.info("Found %d", len(organism.genes))
    LOGGER.info("Loading expression name to gene map")
    name2gene = pyorganism.read_pickle(os.path.join(db_path, "name2gene.pkl"))
    num = sum(1 for gene in name2gene.itervalues() if gene)
    LOGGER.info("Map contains %d names and %d genes (%3.2f%%)", len(name2gene),
            num, 100.0 * num / len(name2gene))
    LOGGER.info("Loading blattner number to gene map")
    blattner2gene = pyorganism.read_pickle(os.path.join(db_path, "blattner2gene.pkl"))
    num = sum(1 for gene in blattner2gene.itervalues() if gene)
    LOGGER.info("Map contains %d blattner numbers and %d genes (%3.2f%%)", len(blattner2gene),
            num, 100.0 * num / len(blattner2gene))
    for (analysis, data_sets, sig_func) in izip(analyses, analysis_sets, significance_functions):
        LOGGER.info("Analysis '%s':", analysis)
        organism.significant[analysis] = dict()
        for name in data_sets:
            LOGGER.info("    Experiment '%s':", name)
            df = organism.activity[name]
            organism.significant[analysis][name] = sig_func(df, name2gene, blattner2gene)

def main(db_path):
    (base_dir, version) = os.path.split(db_path)
    if not version:
        (base_dir, version) = os.path.split(os.path.dirname(db_path))
    ecoli = pyorganism.Organism(name="E. coli K12")
    if CONFIG["continuous"]:
        pass
    else:
        load_discrete(db_path, ecoli, CONFIG["data_paths"], CONFIG["data_load"])
        prepare_discrete(db_path, ecoli, CONFIG["analyses"], CONFIG["data_sets"],
                CONFIG["significance"])
        results = ResultManager(os.path.join(base_dir, "control_results.h5"),
                "/Discrete")
    # general parallel setup using IPython.parallel
    rc = Client()
    dv = rc.direct_view()
#    lv = rc.load_balanced_view()
    lv = None
    dv.execute("import pyorganism;import pyorganism.regulation as pyreg", block=True)
    for (analysis, data, random_num, jackknife_fraction, jackknife_num) in izip(
            CONFIG["analyses"], CONFIG["data_sets"], CONFIG["random_num"], CONFIG["jackknife_fraction"], CONFIG["jackknife_num"]):
        LOGGER.info("*" * 79)
        LOGGER.info(analysis.upper())
        LOGGER.info("*" * 79)
        analysis_function = ANAL_RESOLVE.get(analysis, StandardError("unknown name"))
        ANAL_PREP.get(analysis, StandardError("unknown name"))(db_path, ecoli, dv)
        for name in data:
            analysis_function(analysis, version, name, ecoli, random_num,
                    jackknife_fraction, jackknife_num, results, dv, lv)
    results.finalize()


if __name__ == "__main__":
    main(sys.argv[1])

