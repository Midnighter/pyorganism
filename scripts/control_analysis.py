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
import pandas

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
    return active

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
# other parameters, per analysis
        random_num=[
                1E4,
                1E4
        ],
        jackknife_fraction=[
                0.1,
                0.1
        ],
        jackknife_num=[
                4,
                4
        ]
)

def digital_analysis(version, data_name, trn, active, random_num, jackknife_fraction,
        jackknife_num, results, dv, lv=None):
    LOGGER.info("*" * 79)
    LOGGER.info(data_name)
    LOGGER.info("*" * 79)
    LOGGER.info("Computing digital control")
    dig = pyreg.digital_control(trn, active)
    LOGGER.info("Computing digital ctc")
#    (dig_ctc, sample) = pyreg.digital_ctc(trn, active, random_num, return_sample=True)
    (dig_ctc, sample) = pp.digital_ctc(dv, trn,
        active, random_num, return_sample=True, lb_view=lv)
#        active, random_num, return_sample=True)
    LOGGER.info("Computing digital ctc robustness")
    LOGGER.setLevel(logging.WARN)
    jack_num = int(numpy.floor(len(active) * jackknife_fraction))
    robust = list()
    for i in range(jackknife_num):
        jacked = copy(active)
        jacked = [jacked.pop(numpy.random.randint(len(jacked))) for j in
            range(jack_num)]
#        z_score = pyreg.digital_ctc(trn, active, random_num)
        z_score = pp.digital_ctc(dv, trn, active, random_num)
#        z_score = pp.digital_ctc(dv, trn, active, random_num, lb_view=lv)
        robust.append(z_score)
    LOGGER.setLevel(logging.INFO)
    robust = numpy.array(robust, dtype=float)
    results.append(version, "digital", False, dig, dig_ctc, data_name,
            sample, robust)

ANAL_RESOLVE = dict(
        digital=digital_analysis,
#        analog=analog_analysis,
#        metabolic=metabolic_analysis
)

def load_time_resolved(db_path, organism, experiment_paths):
    pass

def load_discrete(db_path, organism, experiment_paths, loading_funcs):
    LOGGER.info("Loading discrete expression data")
    for ((path, name), load_data) in izip(experiment_paths, loading_funcs):
        LOGGER.info("%s: '%s'", name, path)
        organism.activity[name] = load_data(path)

def prepare_discrete(db_path, organism, analyses):
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
    for data_sets in analyses:
        for name in data_sets:
            df = organism.activity[name]
            organism.significant[name] = simple_discrete(df, name2gene,
                    blattner2gene)
            LOGGER.info("Experiment '%s': %d DEGs", name, len(organism.significant[name]))

def prepare_digital(db_path, organism):
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

def main(db_path):
    (base_dir, version) = os.path.split(db_path)
    if not version:
        (base_dir, version) = os.path.split(os.path.dirname(db_path))
    ecoli = pyorganism.Organism(name="E. coli K12")
    if CONFIG["continuous"]:
        load_time_resolved(db_path, ecoli, CONFIG["data_paths"])
        results = ResultManager(os.path.join(base_dir, "control_results.h5"),
                "/Continuous")
    else:
        load_discrete(db_path, ecoli, CONFIG["data_paths"], CONFIG["data_load"])
        prepare_discrete(db_path, ecoli, CONFIG["data_sets"])
        results = ResultManager(os.path.join(base_dir, "control_results.h5"),
                "/Discrete")
    # general parallel setup using IPython.parallel
    rc = Client()
    dv = rc.direct_view()
    lv = rc.load_balanced_view()
    dv.execute("import pyorganism;import pyorganism.regulation as pyreg", block=True)
    for (analysis, data, random_num, jackknife_fraction, jackknife_num) in izip(
            CONFIG["analyses"], CONFIG["data_sets"], CONFIG["random_num"], CONFIG["jackknife_fraction"], CONFIG["jackknife_num"]):
        analysis_function = ANAL_RESOLVE.get(analysis, StandardError("unknown name"))
        prepare_digital(db_path, ecoli)
        for name in data:
            analysis_function(version, name, ecoli.trn, ecoli.significant[name], random_num,
                    jackknife_fraction, jackknife_num, results, dv, lv)
        break
    results.finalize()


if __name__ == "__main__":
    main(sys.argv[1])

