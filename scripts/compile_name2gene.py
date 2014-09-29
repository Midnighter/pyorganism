#!/usr/bin/env python
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

from pyorganism.io.microarray import read_microarray


LOGGER = logging.getLogger()
LOGGER.addHandler(logging.StreamHandler())
LOGGER.setLevel(logging.INFO)

VERSION = "default"


single_dir = "Expression/LZ41-LZ54_single_knockouts"
double_dir = "Expression/LZ41-LZ54_double_knockouts"
CONFIG = dict(
        continuous=False,
        data_paths=[
            (os.path.join(single_dir, "LZ41-LZ54.tsv"), "wt"),
            (os.path.join(single_dir, "LZ41_d_fis-LZ54_d_fis.tsv"), "fis"),
            (os.path.join(single_dir, "LZ41_d_hns-LZ54_d_hns.tsv"), "hns"),
            (os.path.join(single_dir, "LZ41-LZ41_d_fis.tsv"), "wt-fis-low"),
            (os.path.join(single_dir, "LZ54-LZ54_d_fis.tsv"), "wt-fis-high"),
            (os.path.join(single_dir, "LZ41-LZ41_d_hns.tsv"), "wt-hns-low"),
            (os.path.join(single_dir, "LZ54-LZ54_d_hns.tsv"), "wt-hns-high"),
            (os.path.join(double_dir, "LZ41_d_fis_hns-LZ54_d_fis_hns.tsv"), "fis-hns"),
            (os.path.join(double_dir, "LZ41-LZ41_d_fis_hns.tsv"), "wt-fis-hns-low"),
            (os.path.join(double_dir, "LZ54-LZ54_d_fis_hns.tsv"),
                    "wt-fis-hns-high"),
        ],
        data_load=[
            read_microarray,
            read_microarray,
            read_microarray,
            read_microarray,
            read_microarray,
            read_microarray,
            read_microarray,
            read_microarray,
            read_microarray,
            read_microarray
        ],
        threshold=80
)


def compile_names(experiment_paths, loading_funcs):
    LOGGER.info("Loading discrete expression data")
    data_frames = [load_data(path)\
            for ((path, name), load_data) in izip(experiment_paths, loading_funcs)]
    return pandas.concat(data_frames, ignore_index=True)

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
            targets.append(gene.bnumber)
            indeces.append(i)
        for name in gene.synonyms:
            targets.append(name)
            indeces.append(i)
#        if gene.product:
#            targets.append(gene.product.name)
#            indeces.append(i)
#            for name in gene.product.synonyms:
#                targets.append(name)
#                indeces.append(i)
#        if gene.regulatory_product:
#            targets.append(gene.regulatory_product.name)
#            indeces.append(i)
#            for name in gene.regulatory_product.synonyms:
#                targets.append(name)
#                indeces.append(i)
    return pyorganism.FindObject(genes, targets=targets, indeces=indeces)

def manual_name_updates(name2gene):
# [b1500](http://regulondb.ccg.unam.mx/gene?term=ECK120003379&organism=ECK12&format=jsp&type=gene)
    gene = pyreg.Gene.get("ECK120003379", None, VERSION)
    name2gene["b1500"] = gene
# [b0609](http://regulondb.ccg.unam.mx/gene?term=ECK120002943&organism=ECK12&format=jsp&type=gene)
    gene = pyreg.Gene.get("ECK120002943", None, VERSION)
    name2gene["b0609"] = gene
# [b0332](http://www.ncbi.nlm.nih.gov/gene/?term=944992) is not recorded in RegulonDB.
# [b3975](http://regulondb.ccg.unam.mx/gene?term=G7818&type=gene&format=jsp)
    gene = pyreg.Gene.get("ECK120004318", None, VERSION)
    name2gene["b3975"] = gene
# [b2084](http://regulondb.ccg.unam.mx/gene?term=G7121&type=gene&format=jsp)
    gene = pyreg.Gene.get("ECK120003706", None, VERSION)
    name2gene["b2084"] = gene
# [b0322](http://regulondb.ccg.unam.mx/gene?term=G0-10550&type=gene&format=jsp)
    gene = pyreg.Gene.get("ECK120026444", None, VERSION)
    name2gene["b0322"] = gene
# [b0309](http://regulondb.ccg.unam.mx/gene?term=ECK120002798&organism=ECK12&format=jsp&type=gene)
    gene = pyreg.Gene.get("ECK120002798", None, VERSION)
    name2gene["b0309"] = gene
# [b2596](http://regulondb.ccg.unam.mx/gene?term=G7353&type=gene&format=jsp)
    gene = pyreg.Gene.get("ECK120003922", None, VERSION)
    name2gene["b2596"] = gene
# [b1364](http://regulondb.ccg.unam.mx/gene?term=ECK120003282&organism=ECK12&format=jsp&type=gene)
    gene = pyreg.Gene.get("ECK120003282", None, VERSION)
    name2gene["b1364"] = gene
# [b4091](http://porteco.org/AjaxSearch.jsp?searchString=b4091) is not recorded in RegulonDB.
    gene = pyreg.Gene.get("ECK120003116", None, VERSION)
    name2gene["ycdF"] = gene
# [ycdF](http://regulondb.ccg.unam.mx/gene?term=ECK120003116&organism=ECK12&format=jsp&type=gene)
    gene = pyreg.Gene.get("ECK120003116", None, VERSION)
    name2gene["ycdF"] = gene
# [b1903](http://regulondb.ccg.unam.mx/gene?term=G7034&type=gene&format=jsp)
    gene = pyreg.Gene.get("ECK120003620", None, VERSION)
    name2gene["b1903"] = gene
# [b0501](http://regulondb.ccg.unam.mx/gene?term=G6272&type=gene&format=jsp)
    gene = pyreg.Gene.get("ECK120002880", None, VERSION)
    name2gene["b0501"] = gene
# [b0165](http://regulondb.ccg.unam.mx/gene?term=ECK120002711&organism=ECK12&format=jsp&type=gene)
    gene = pyreg.Gene.get("ECK120002711", None, VERSION)
    name2gene["b0165"] = gene
# [rhsC_1](http://regulondb.ccg.unam.mx/gene?term=ECK120000839&organism=ECK12&format=jsp&type=gene)
    gene = pyreg.Gene.get("ECK120000839", None, VERSION)
    name2gene["rhsC_1"] = gene
# [rhsC_2](http://regulondb.ccg.unam.mx/gene?term=ECK120000839&organism=ECK12&format=jsp&type=gene)
    gene = pyreg.Gene.get("ECK120000839", None, VERSION)
    name2gene["rhsC_2"] = gene
# [b1354](http://regulondb.ccg.unam.mx/gene?term=G6678&type=gene&format=jsp)
    gene = pyreg.Gene.get("ECK120003273", None, VERSION)
    name2gene["b1354"] = gene
# [b2651](http://regulondb.ccg.unam.mx/gene?term=ECK120003954&organism=ECK12&format=jsp&type=gene)
    gene = pyreg.Gene.get("ECK120003954", None, VERSION)
    name2gene["b2651"] = gene
# [b1052](http://regulondb.ccg.unam.mx/gene?term=ECK120003149&organism=ECK12&format=jsp&type=gene)
    gene = pyreg.Gene.get("ECK120003149", None, VERSION)
    name2gene["b1052"] = gene
# [b0395](http://www.ncbi.nlm.nih.gov/gene/?term=949074) is not recorded in RegulonDB.
# [b0725](http://regulondb.ccg.unam.mx/gene?term=ECK120002988&organism=ECK12&format=jsp&type=gene)
    gene = pyreg.Gene.get("ECK120002988", None, VERSION)
    name2gene["b0725"] = gene
# [b3837](http://porteco.org/AjaxSearch.jsp?searchString=b3837) is not recorded in RegulonDB.
# [b3007](http://regulondb.ccg.unam.mx/gene?term=G7562&type=gene&format=jsp)
    gene = pyreg.Gene.get("ECK120004121", None, VERSION)
    name2gene["b3007"] = gene
# [b3004](http://regulondb.ccg.unam.mx/gene?term=G7561&type=gene&format=jsp)
    gene = pyreg.Gene.get("ECK120004120", None, VERSION)
    name2gene["b3004"] = gene

def verify_union(name2gene, blattner2gene, full_data):
    conflicts = 0
    matched = 0
    pairs = set(full_data[["name", "blattner"]].itertuples(index=False))
    for (name, blattner) in pairs:
        name_gene = name2gene.get(name)
        blattner_gene = blattner2gene.get(blattner)
        if name_gene or blattner_gene:
            matched += 1
        if blattner_gene and name_gene and name_gene.unique_id != blattner_gene.unique_id:
            LOGGER.warn("conflict for %s: %s and %s: %s", name,
                    name_gene.unique_id, blattner, blattner_gene.unique_id)
            conflicts += 1
    LOGGER.info("Using the union of the maps %d genes were found (%3.2f%%)"\
            " involving %d conflicts", matched,
            100.0 * matched / len(pairs), conflicts)


def compile_name2gene(objects_path):
    LOGGER.info("{0:*^78s}".format("Compile Gene Name Map"))
    full_data = compile_names(CONFIG["data_paths"], CONFIG["data_load"])
    version = os.path.basename(objects_path)
    if not version:
        version = os.path.basename(os.path.dirname(objects_path))
    global VERSION
    VERSION = version
    LOGGER.info("{0:*^78s}".format(version))
    LOGGER.info("Loading genes")
    genes = pyorganism.read_pickle(os.path.join(objects_path, "genes.pkl"))
    LOGGER.info("Finding gene names")
    names = set(full_data["name"].unique())
    if numpy.nan in names:
        names.remove(numpy.nan)
    names = sorted(names)
    finder = pyorganism.FindObject(genes, "name")
    (name2gene, name_gap) = compile_feature_map(genes, names, finder)
    synonyms = synonym_finder(genes)
    LOGGER.info("Finding gene names by synonyms")
    name_gap = extend_feature_map(name2gene, name_gap, synonyms)
    LOGGER.info("Missing %d gene names", len(name_gap))
    LOGGER.info("Fuzzy search of gene names (threshold %d%%)", CONFIG["threshold"])
    name_gap = fuzzy_extension(name2gene, name_gap, finder, CONFIG["threshold"])
#    LOGGER.info("Fuzzy search of gene names by synonyms (threshold %d%%)", CONFIG["threshold"])
#    name_gap = fuzzy_extension(name2gene, name_gap, synonyms, CONFIG["threshold"])
    manual_name_updates(name2gene)
    num = sum(1 for gene in name2gene.itervalues() if gene)
    LOGGER.info("Final map contains %d names and %d genes (%3.2f%%)", len(name2gene),
            num, 100.0 * num / len(name2gene))
    LOGGER.info("Finding gene blattner numbers")
    bnumbers = set(full_data["blattner"].unique())
    if numpy.nan in bnumbers:
        bnumbers.remove(numpy.nan)
    bnumbers = sorted(bnumbers)
    gene_finder = pyorganism.FindObject(genes, "bnumber")
    (blattner2gene, blattner_gap) = compile_feature_map(genes, bnumbers, gene_finder)
    LOGGER.info("Finding gene blattner numbers by synonyms")
    blattner_gap = extend_feature_map(blattner2gene, blattner_gap, synonyms)
    LOGGER.info("Missing %d gene blattner numbers", len(blattner_gap))
    LOGGER.info("Fuzzy search of gene blattner numbers (threshold %d%%)", CONFIG["threshold"])
    blattner_gap = fuzzy_extension(blattner2gene, blattner_gap, finder, CONFIG["threshold"])
#    LOGGER.info("Fuzzy search of gene blattner numbers by synonyms (threshold %d%%)", CONFIG["threshold"])
#    blattner_gap = fuzzy_extension(blattner2gene, blattner_gap, synonyms, CONFIG["threshold"])
    num = sum(1 for gene in blattner2gene.itervalues() if gene)
    LOGGER.info("Final map contains %d blattner numbers and %d genes (%3.2f%%)",
            len(blattner2gene), num, 100.0 * num / len(blattner2gene))
    verify_union(name2gene, blattner2gene, full_data)
    pyorganism.write_pickle(name2gene, os.path.join(objects_path, "name2gene.pkl"))
    pyorganism.write_pickle(blattner2gene, os.path.join(objects_path, "blattner2gene.pkl"))

if __name__ == "__main__":
    if len(sys.argv) != 2:
        LOGGER.critical("%s <RegulonDB objects path>", sys.argv[0])
        sys.exit(2)
    else:
        compile_name2gene(sys.argv[1])

