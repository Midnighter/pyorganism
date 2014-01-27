# -*- coding: utf-8 -*-


import sys
import os
import logging

import numpy
import pandas

import pyorganism
import pyorganism.regulation as pyreg

from itertools import izip

from fuzzywuzzy import fuzz

from pyorganism.io.microarray import read_interpolated


LOGGER = logging.getLogger()
LOGGER.addHandler(logging.StreamHandler())
LOGGER.setLevel(logging.INFO)


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
        threshold=80
)


def compile_names(experiment_paths, loading_funcs):
    LOGGER.info("Loading interpolated RNAseq data")
    data_frames = [load_data(path)\
            for ((path, name), load_data) in izip(experiment_paths, loading_funcs)]
    names = set()
    for df in data_frames:
        names.update(str(x) for x in df.index)
    # no need to test for missing values, index must be complete
    return sorted(names)

def compile_feature_map(genes, features, gene_finder):
    feature2gene = dict()
    gap = list()
    for name in features:
        try:
            feature2gene[name] = gene_finder(name)
        except IndexError:
            gap.append(name)
    LOGGER.info("Found %d/%d (%3.2f%%)", len(feature2gene), len(features),
            100.0 * len(feature2gene) / len(features))
    return (feature2gene, gap)

def extend_feature_map(feature2gene, gap, gene_finder):
    found = list()
    for name in gap:
        try:
            feature2gene[name] = gene_finder(name)
            found.append(name)
        except IndexError:
            continue
    gap = set(gap).difference(set(found))
    LOGGER.info("Found %d additional genes", len(found))
    return gap

def fuzzy_extension(feature2gene, gap, gene_finder, thresh=80):
    found = list()
    for name in gap:
        try:
            (gene, match, score) = gene_finder.fuzzy_search(name,
                    threshold=thresh, scorer=fuzz.QRatio)
            feature2gene[name] = gene
            found.append(name)
        except IndexError:
            LOGGER.debug("'%s' not found", name)
            feature2gene[name] = None
    gap = set(gap).difference(set(found))
    LOGGER.info("Found %d additional genes", len(found))
    return gap

def synonym_finder(genes):
    targets = list()
    indeces = list()
    for (i, gene) in enumerate(genes):
        if gene.bnumber:
            targets.append(gene.bnumber.lower())
            indeces.append(i)
        for name in gene.synonyms:
            targets.append(name.lower())
            indeces.append(i)
        if gene.product:
            targets.append(gene.product.name.lower())
            indeces.append(i)
            for name in gene.product.synonyms:
                targets.append(name.lower())
                indeces.append(i)
        if gene.regulatory_product:
            targets.append(gene.regulatory_product.name.lower())
            indeces.append(i)
            for name in gene.regulatory_product.synonyms:
                targets.append(name.lower())
                indeces.append(i)
    return pyorganism.FindObject(genes, targets=targets, indeces=indeces)

def manual_name_updates(feature2gene):
# [waak](http://www.ecogene.org/old/geneinfo.php?eg_id=EG11423)
    gene = pyreg.Gene.get("ECK120001387")
    feature2gene["waak"] = gene
# [ybbv](http://ecoliwiki.net/colipedia/index.php/ybbV:Quickview)
    gene = pyreg.Gene.get("ECK120002943")
    feature2gene["ybbv"] = gene

def compile_feature2gene(objects_path):
    version = os.path.basename(objects_path)
    if not version:
        version = os.path.basename(os.path.dirname(objects_path))
    names = compile_names(CONFIG["data_paths"], CONFIG["data_load"])
    LOGGER.info("Loading genes")
    genes = pyorganism.read_pickle(os.path.join(objects_path, "genes.pkl"))
    LOGGER.info("Finding gene names")
    finder = pyorganism.FindObject(genes,
            targets=[gene.name.lower() for gene in genes])
    (feature2gene, gap) = compile_feature_map(genes, names, finder)
    synonyms = synonym_finder(genes)
    LOGGER.info("Finding gene names by synonyms")
    gap = extend_feature_map(feature2gene, gap, synonyms)
    LOGGER.info("Missing %d gene names", len(gap))
    LOGGER.info("Fuzzy search of gene names (threshold %d%%)", CONFIG["threshold"])
    LOGGER.setLevel(logging.DEBUG)
    gap = fuzzy_extension(feature2gene, gap, finder, CONFIG["threshold"])
    LOGGER.info("Fuzzy search of gene names by synonyms (threshold %d%%)", CONFIG["threshold"])
    gap = fuzzy_extension(feature2gene, gap, synonyms, CONFIG["threshold"])
    manual_name_updates(feature2gene)
    num = sum(1 for gene in feature2gene.itervalues() if gene)
    LOGGER.info("Final map contains %d names and %d genes (%3.2f%%)", len(feature2gene),
            num, 100.0 * num / len(feature2gene))
    pyorganism.write_pickle(feature2gene, os.path.join(objects_path, "feature2gene.pkl"))

if __name__ == "__main__":
    if len(sys.argv) != 2:
        LOGGER.critical("%s <RegulonDB objects path>", sys.argv[0])
        sys.exit(2)
    else:
        compile_feature2gene(sys.argv[1])

