# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Creating a mapping from gene expression data to their respective genes

# <codecell>

import os

import numpy
import pandas

import pyorganism
import pyorganism.regulation as pyreg

from glob import glob

# <markdowncell>

# We extract all expression data identifiers from the microarray data.

# <codecell>

ecoli = pyorganism.Organism(name="E. coli K12")

# <codecell>

base_dir = "Data/Expression/LZ41-LZ54_single_knockouts"

# <codecell>

ecoli.read_activity(os.path.join(base_dir, "LZ41-LZ54.tsv"), "wt");

# <codecell>

ecoli.read_activity(os.path.join(base_dir, "LZ41_d_fis-LZ54_d_fis.tsv"), "fis");

# <codecell>

ecoli.read_activity(os.path.join(base_dir, "LZ41_d_hns-LZ54_d_hns.tsv"), "hns");

# <codecell>

ecoli.read_activity(os.path.join(base_dir, "LZ41-LZ41_d_fis.tsv"), "wt-fis-low");

# <codecell>

ecoli.read_activity(os.path.join(base_dir, "LZ54-LZ54_d_fis.tsv"), "wt-fis-high");

# <codecell>

ecoli.read_activity(os.path.join(base_dir, "LZ41-LZ41_d_hns.tsv"), "wt-hns-low");

# <codecell>

ecoli.read_activity(os.path.join(base_dir, "LZ54-LZ54_d_hns.tsv"), "wt-hns-high");

# <codecell>

base_dir = "Data/Expression/LZ41-LZ54_double_knockouts"

# <codecell>

ecoli.read_activity(os.path.join(base_dir, "LZ41_d_fis_hns-LZ54_d_fis_hns.tsv"), "fis-hns");

# <codecell>

ecoli.read_activity(os.path.join(base_dir, "LZ41-LZ41_d_fis_hns.tsv"), "wt-fis-hns-low");

# <codecell>

ecoli.read_activity(os.path.join(base_dir, "LZ54-LZ54_d_fis_hns.tsv"), "wt-fis-hns-high");

# <markdowncell>

# How many unique entities are in the expression data?

# <codecell>

expression = set()
for data in ecoli.activity.itervalues():
    expression.update([x[0] for x in data])
expression = sorted(expression)
print len(expression)

# <markdowncell>

# Select a database version.

# <codecell>

versions = sorted(glob("Data/RegulonDBObjects/[0-9].[0-9]"))
for (i, ver) in enumerate(versions):
    print "{0:2d} -> {1}".format(i, ver)

# <codecell>

index = 16
path_ver = versions[index]
num_ver = os.path.basename(path_ver)
print num_ver
print path_ver

# <markdowncell>

# Load all genes.

# <codecell>

genes = pyorganism.read_pickle(os.path.join(path_ver, "genes.pkl"))
print len(genes)

# <codecell>

gene_finder = pyorganism.FindObject(genes, "name")

# <codecell>

name2gene = dict()
gap = list()
for name in expression:
    try:
        name2gene[name] = gene_finder(name)
    except IndexError as err:
        gap.append(name)
print "found", len(name2gene), "/", len(expression)
print "missing", len(gap)

# <codecell>

targets = list()
indeces = list()
for (i, gene) in enumerate(genes):
    targets.append(gene.name)
    indeces.append(i)
    if gene.bnumber:
        targets.append(gene.bnumber)
        indeces.append(i)
    for name in gene.synonyms:
        targets.append(name)
        indeces.append(i)
    if gene.product:
        targets.append(gene.product.name)
        indeces.append(i)
        for name in gene.product.synonyms:
            targets.append(name)
            indeces.append(i)
    if gene.regulatory_product:
        targets.append(gene.regulatory_product.name)
        indeces.append(i)
        for name in gene.regulatory_product.synonyms:
            targets.append(name)
            indeces.append(i)

# <codecell>

gene_finder = pyorganism.FindObject(genes, targets=targets, indeces=indeces)

# <codecell>

found = list()
for name in gap:
    try:
        name2gene[name] = gene_finder(name)
        found.append(name)
    except IndexError as err:
        continue
gap = set(gap).difference(set(found))
print "found", len(name2gene), "/", len(expression)
print "missing", len(gap)

# <codecell>

from fuzzywuzzy import fuzz

# <markdowncell>

# From all expression data, which are the entries without association?

# <codecell>

for name in gap:
    try:
        (gene, match, score) = gene_finder.fuzzy_search(name, threshold=80, scorer=fuzz.QRatio)
        print name, match, score
        name2gene[name] = gene
    except IndexError as err:
        name2gene[name] = None
        print err

# <markdowncell>

# Unknown names (version dependent):

# <markdowncell>

# [b1500](http://regulondb.ccg.unam.mx/gene?term=ECK120003379&organism=ECK12&format=jsp&type=gene)

# <codecell>

gene = pyreg.Gene.get("ECK120003379", None)
name2gene["b1500"] = gene
print gene

# <markdowncell>

# [b0609](http://regulondb.ccg.unam.mx/gene?term=ECK120002943&organism=ECK12&format=jsp&type=gene)

# <codecell>

gene = pyreg.Gene.get("ECK120002943", None)
name2gene["b0609"] = gene
print gene

# <markdowncell>

# [b0332](http://www.ncbi.nlm.nih.gov/gene/?term=944992) is not recorded in RegulonDB.

# <markdowncell>

# [b3975](http://regulondb.ccg.unam.mx/gene?term=G7818&type=gene&format=jsp)

# <codecell>

gene = pyreg.Gene.get("ECK120004318", None)
name2gene["b3975"] = gene
print gene

# <markdowncell>

# [b2084](http://regulondb.ccg.unam.mx/gene?term=G7121&type=gene&format=jsp)

# <codecell>

gene = pyreg.Gene.get("ECK120003706", None)
name2gene["b2084"] = gene
print gene

# <markdowncell>

# [b0322](http://regulondb.ccg.unam.mx/gene?term=G0-10550&type=gene&format=jsp)

# <codecell>

gene = pyreg.Gene.get("ECK120026444", None)
name2gene["b0322"] = gene
print gene

# <markdowncell>

# [b2391](http://www.ncbi.nlm.nih.gov/gene/?term=946863) is not recorded in RegulonDB.

# <markdowncell>

# [b0309](http://regulondb.ccg.unam.mx/gene?term=ECK120002798&organism=ECK12&format=jsp&type=gene)

# <codecell>

gene = pyreg.Gene.get("ECK120002798", None)
name2gene["b0309"] = gene
print gene

# <markdowncell>

# [b2596](http://regulondb.ccg.unam.mx/gene?term=G7353&type=gene&format=jsp)

# <codecell>

gene = pyreg.Gene.get("ECK120003922", None)
name2gene["b2596"] = gene
print gene

# <markdowncell>

# [b1364](http://regulondb.ccg.unam.mx/gene?term=ECK120003282&organism=ECK12&format=jsp&type=gene)

# <codecell>

gene = pyreg.Gene.get("ECK120003282", None)
name2gene["b1364"] = gene
print gene

# <markdowncell>

# [b4091](http://porteco.org/AjaxSearch.jsp?searchString=b4091) is not recorded in RegulonDB.

# <markdowncell>

# [ycdF](http://regulondb.ccg.unam.mx/gene?term=ECK120003116&organism=ECK12&format=jsp&type=gene)

# <codecell>

gene = pyreg.Gene.get("ECK120003116", None)
name2gene["ycdF"] = gene
print gene

# <markdowncell>

# [b1903](http://regulondb.ccg.unam.mx/gene?term=G7034&type=gene&format=jsp)

# <codecell>

gene = pyreg.Gene.get("ECK120003620", None)
name2gene["b1903"] = gene
print gene

# <markdowncell>

# [b0501](http://regulondb.ccg.unam.mx/gene?term=G6272&type=gene&format=jsp)

# <codecell>

gene = pyreg.Gene.get("ECK120002880", None)
name2gene["b0501"] = gene
print gene

# <markdowncell>

# [b0165](http://regulondb.ccg.unam.mx/gene?term=ECK120002711&organism=ECK12&format=jsp&type=gene)

# <codecell>

gene = pyreg.Gene.get("ECK120002711", None)
name2gene["b0165"] = gene
print gene

# <markdowncell>

# [rhsC_1](http://regulondb.ccg.unam.mx/gene?term=ECK120000839&organism=ECK12&format=jsp&type=gene)

# <codecell>

gene = pyreg.Gene.get("ECK120000839", None)
name2gene["rhsC_1"] = gene
print gene

# <markdowncell>

# [rhsC_2](http://regulondb.ccg.unam.mx/gene?term=ECK120000839&organism=ECK12&format=jsp&type=gene)

# <codecell>

gene = pyreg.Gene.get("ECK120000839", None)
name2gene["rhsC_2"] = gene
print gene

# <markdowncell>

# [b1354](http://regulondb.ccg.unam.mx/gene?term=G6678&type=gene&format=jsp)

# <codecell>

gene = pyreg.Gene.get("ECK120003273", None)
name2gene["b1354"] = gene
print gene

# <markdowncell>

# [b2651](http://regulondb.ccg.unam.mx/gene?term=ECK120003954&organism=ECK12&format=jsp&type=gene)

# <codecell>

gene = pyreg.Gene.get("ECK120003954", None)
name2gene["b2651"] = gene
print gene

# <markdowncell>

# [b1052](http://regulondb.ccg.unam.mx/gene?term=ECK120003149&organism=ECK12&format=jsp&type=gene)

# <codecell>

gene = pyreg.Gene.get("ECK120003149", None)
name2gene["b1052"] = gene
print gene

# <markdowncell>

# [b0395](http://www.ncbi.nlm.nih.gov/gene/?term=949074) is not recorded in RegulonDB.

# <markdowncell>

# [b0725](http://regulondb.ccg.unam.mx/gene?term=ECK120002988&organism=ECK12&format=jsp&type=gene)

# <codecell>

gene = pyreg.Gene.get("ECK120002988", None)
name2gene["b0725"] = gene
print gene

# <markdowncell>

# [b3837](http://porteco.org/AjaxSearch.jsp?searchString=b3837) is not recorded in RegulonDB.

# <markdowncell>

# [b3007](http://regulondb.ccg.unam.mx/gene?term=G7562&type=gene&format=jsp)

# <codecell>

gene = pyreg.Gene.get("ECK120004121", None)
name2gene["b3007"] = gene
print gene

# <markdowncell>

# [b3004](http://regulondb.ccg.unam.mx/gene?term=G7561&type=gene&format=jsp)

# <codecell>

gene = pyreg.Gene.get("ECK120004120", None)
name2gene["b3004"] = gene
print gene

# <headingcell level=2>

# Storing and exchange of information between versions

# <codecell>

eck12 = list()
for name in expression:
    gene = name2gene.get(name)
    eck12.append(gene.unique_id if gene else "")
df = pandas.DataFrame(data=eck12, index=expression, columns=(num_ver,))
df.sort(inplace=True)

# <codecell>

set(type(item) for item in df)

# <codecell>

(df == "").sum()

# <codecell>

filename = os.path.join("Data/RegulonDBObjects/unified_name2gene.h5")
if os.path.exists(filename):
    mapping = pandas.read_hdf(filename, "/mapping")
    mapping[num_ver] = pandas.Series(df[num_ver], index=mapping.index)
else:
    df.to_hdf(filename, "/mapping", format="table")
    mapping = pandas.read_hdf(filename, "/mapping")

# <codecell>

set(type(item) for item in mapping)

# <markdowncell>

# Removed for now, let's see the differences between versions.
# 
#     columns = list(mapping.columns)
#     columns.remove(num_ver)
#     for col in columns:
#         # this and this_mask may change during loop
#         this = mapping[num_ver]
#         this_mask = (this == "")
#         other = mapping[col]
#         other_mask = (other == "")
#         for (index, value) in this_mask.iteritems():
#             if value:
#                 if not other_mask.loc[index]:
#                     unique_id = other.loc[index]
#                     if unique_id in pyreg.Gene:
#                         this.loc[index] = unique_id

# <codecell>

for col in mapping.columns:
    print "Version", col, "Unmapped names:", (mapping[col] == "").sum()

# <codecell>

faulty_indeces = [row[0] for row in mapping.itertuples() if not all(x == row[1] for x in row[1:])]
mapping.loc[faulty_indeces]

# <codecell>

mapping.to_hdf(filename, "/mapping", format="table")

