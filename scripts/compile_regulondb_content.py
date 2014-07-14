#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
=========================
Compile RegulonDB Content
=========================

:Authors:
    Moritz Emanuel Beber
:Date:
    2013-04-02
:Copyright:
    Copyright |c| 2013 Jacobs University Bremen, all rights reserved.
:File:
    compile_regulondb_content.py

.. |c| unicode:: U+A9
"""


import logging
LOGGER = logging.getLogger()
LOGGER.addHandler(logging.StreamHandler())
LOGGER.setLevel(logging.INFO)


import sys
import os

import pyorganism
import pyorganism.regulation as pyreg
import pyorganism.io.regulondb as regdb


def compile_genes(base_dir):
    LOGGER.info("Compiling genes...")
    genes = regdb.read_genes(os.path.join(base_dir, "gene.xml"))
    norm = len(genes)
    LOGGER.info("Found {0:d}".format(len(genes)))
    if norm == 0:
        LOGGER.error("Failed to compile genes.")
        return genes
    LOGGER.info("Compiling additional information...")
    regdb.update_gene_synonyms(os.path.join(base_dir, "object_synonym.xml"))
    regdb.update_gene_external(os.path.join(base_dir, "object_external_db_link.xml"))
    norm = float(norm)
    num = sum(1 for gene in genes if gene.name)
    LOGGER.info("Found {0:d} genes with names ({1:.2%})".format(num, num / norm))
    num = sum(1 for gene in genes if gene.bnumber)
    LOGGER.info("Found {0:d} genes with Blattner number ({1:.2%})".format(num, num / norm))
    num = sum(1 for gene in genes if gene.synonyms)
    LOGGER.info("Found {0:d} genes with synonyms ({1:.2%})".format(num, num / norm))
    return genes

def compile_products(base_dir):
    LOGGER.info("Compiling products...")
    products = regdb.read_products(os.path.join(base_dir, "product.xml"))
    norm = len(products)
    LOGGER.info("Found {0:d}".format(norm))
    if norm == 0:
        LOGGER.error("Failed to compile products.")
        return products
    LOGGER.info("Compiling additional information...")
    regdb.update_synonyms(os.path.join(base_dir, "object_synonym.xml"),
            pyreg.Product)
    regdb.update_external(os.path.join(base_dir, "object_external_db_link.xml"),
            pyreg.Product)
    norm = float(norm)
    num = sum(1 for prod in products if prod.name)
    LOGGER.info("Found {0:d} products with names ({1:.2%})".format(num, num / norm))
    num = sum(1 for prod in products if prod.synonyms)
    LOGGER.info("Found {0:d} products with synonyms ({1:.2%})".format(num, num / norm))
    num = sum(1 for prod in products if prod.go)
    LOGGER.info("Found {0:d} products with GO ID ({1:.2%})".format(num, num / norm))
    return products

def compile_transcription_factors(base_dir, version):
    # use transcription_factor.xml as it is more complete and only missing from 6.0
    if version == (6, 0):
        LOGGER.error("no transcription_factor.xml file for this version")
        return
    LOGGER.info("Compiling transcription factors...")
    t_factors = regdb.read_transcription_factors( os.path.join(base_dir,
            "transcription_factor.xml"))
    norm = len(t_factors)
    LOGGER.info("Found {0:d}".format(norm))
    if norm == 0:
        LOGGER.error("Failed to compile transcription factors.")
        return t_factors
    regdb.update_synonyms(os.path.join(base_dir, "object_synonym.xml"),
            pyreg.TranscriptionFactor)
    regdb.update_external(os.path.join(base_dir, "object_external_db_link.xml"),
            pyreg.TranscriptionFactor)
    norm = float(norm)
    num = sum(1 for tf in t_factors if tf.name)
    LOGGER.info("Found {0:d} transcription factors with a name ({1:.2%})".format(num, num / norm))
    num = sum(1 for tf in t_factors if tf.synonyms)
    LOGGER.info("Found {0:d} transcription factors with synonyms ({1:.2%})".format(num, num / norm))
    regdb.update_product_transcription_factor_link(os.path.join(base_dir,
            "product_tf_link.xml"))
    return t_factors

def compile_sigma_factors(base_dir):
    LOGGER.info("Compiling sigma factors...")
    sigma_factors = regdb.read_sigma_factors(os.path.join(base_dir, "sigma_tmp.xml"))
    norm = len(sigma_factors)
    LOGGER.info("Found {0:d}".format(norm))
    if norm == 0:
        LOGGER.error("Failed to compile sigma factors.")
        return sigma_factors
    return sigma_factors

def compile_conformations(base_dir):
    LOGGER.info("Compiling conformations...")
    conformations = regdb.read_conformations(os.path.join(base_dir, "conformation.xml"))
    norm = len(conformations)
    LOGGER.info("Found {0:d}".format(norm))
    if norm == 0:
        LOGGER.error("Failed to compile conformations.")
        return conformations
    norm = float(norm)
    num = sum(1 for conf in conformations if conf.final_state)
    LOGGER.info("Found {0:d} conformations with final state ({1:.2%})".format(num, num / norm))
    num = sum(1 for conf in conformations if conf.type)
    LOGGER.info("Found {0:d} conformations with conformation type ({1:.2%})".format(num, num / norm))
    num = sum(1 for conf in conformations if conf.apo_holo == u"Apo")
    LOGGER.info("Found {0:d} apo conformations ({1:.2%})".format(num, num / norm))
    num = sum(1 for conf in conformations if conf.apo_holo == u"Holo")
    LOGGER.info("Found {0:d} holo conformations ({1:.2%})".format(num, num / norm))
    return conformations

def compile_promoters(base_dir):
    LOGGER.info("Compiling promoters...")
    promoters = regdb.read_promoters(os.path.join(base_dir, "promoter.xml"))
    norm = len(promoters)
    LOGGER.info("Found {0:d}".format(norm))
    if norm == 0:
        LOGGER.error("Failed to compile promoters.")
        return promoters
    return promoters

def compile_operons(base_dir):
    LOGGER.info("Compiling operons...")
    operons = regdb.read_operons(os.path.join(base_dir, "operon.xml"))
    norm = len(operons)
    LOGGER.info("Found {0:d}".format(norm))
    if norm == 0:
        LOGGER.error("Failed to compile operons.")
        return operons
    return operons

def compile_transcription_units(base_dir):
    LOGGER.info("Compiling transcription units...")
    t_units = regdb.read_transcription_units(os.path.join(base_dir,
            "transcription_unit.xml"))
    norm = len(t_units)
    LOGGER.info("Found {0:d}".format(norm))
    if norm == 0:
        LOGGER.error("Failed to compile transcription units.")
        return t_units
    regdb.link_tu_genes(os.path.join(base_dir, "tu_gene_link.xml"))
    norm = float(norm)
    num = sum(1 for tu in t_units if tu.genes)
    LOGGER.info("Found {0:d} transcription units that contain genes ({1:.2%})".format(num, num / norm))
    num = sum(1 for tu in t_units if tu.promoter)
    LOGGER.info("Found {0:d} transcription units with associated promoters ({1:.2%})".format(num, num / norm))
    return t_units

def compile_regulation(base_dir):
    LOGGER.info("Compiling regulatory interactions...")
    interactions = regdb.read_regulatory_interactions(os.path.join(base_dir, "regulatory_interaction.xml"))
    norm = len(interactions)
    LOGGER.info("Found {0:d}".format(norm))
    if norm == 0:
        LOGGER.error("Failed to compile regulatory interactions.")
        return interactions
    return interactions


def main(argv):
    if not os.path.exists(argv[0]):
        sys.exit(1)
    if not os.path.exists(argv[1]):
        os.makedirs(argv[1])
    version = argv[2] if len(argv) == 3 else os.path.basename(argv[0])
    if not version:
        version = os.path.basename(os.path.dirname(argv[0]))
    LOGGER.info("{0:*^78s}".format(version))
    version = tuple([int(num) for num in version.split(".")])
    # compile genes and gene products
    gene_file = os.path.join(argv[1], "genes.pkl")
    if os.path.exists(gene_file):
        genes = pyorganism.read_pickle(gene_file)
    else:
        genes = compile_genes(argv[0])
    product_file = os.path.join(argv[1], "products.pkl")
    if os.path.exists(product_file):
        products = pyorganism.read_pickle(product_file)
    else:
        products = compile_products(argv[0])
        regdb.link_gene_product(os.path.join(argv[0], "gene_product_link.xml"))
    tf_file = os.path.join(argv[1], "transcription_factors.pkl")
    if os.path.exists(tf_file):
        t_factors = pyorganism.read_pickle(tf_file)
    else:
        t_factors = compile_transcription_factors(argv[0], version)
        num = sum(1 for gene in genes if isinstance(gene.regulatory_product, pyreg.TranscriptionFactor))
        norm = float(len(genes))
        LOGGER.info("Found {0:d} genes that code for transcription factors({1:.2%})".format(num, num / norm))
    if version >= (7, 2):
        sigma_file = os.path.join(argv[1], "sigma_factors.pkl")
        if os.path.exists(sigma_file):
            sigma_factors = pyorganism.read_pickle(sigma_file)
        else:
            sigma_factors = compile_sigma_factors(argv[0])
    conf_file = os.path.join(argv[1], "conformations.pkl")
    if os.path.exists(conf_file):
        conformations = pyorganism.read_pickle(conf_file)
    else:
        conformations = compile_conformations(argv[0])
    prom_file = os.path.join(argv[1], "promoters.pkl")
    if os.path.exists(prom_file):
        promoters = pyorganism.read_pickle(prom_file)
    else:
        promoters = compile_promoters(argv[0])
    op_file = os.path.join(argv[1], "operons.pkl")
    if os.path.exists(op_file):
        operons = pyorganism.read_pickle(op_file)
    else:
        operons = compile_operons(argv[0])
        regdb.update_operons(operons, promoters, genes)
    tu_file = os.path.join(argv[1], "transcription_units.pkl")
    if os.path.exists(tu_file):
        t_units = pyorganism.read_pickle(tu_file)
    else:
        t_units = compile_transcription_units(argv[0])
    inter_file = os.path.join(argv[1], "interactions.pkl")
    if os.path.exists(inter_file):
        interactions = pyorganism.read_pickle(inter_file)
    else:
        interactions = compile_regulation(argv[0])
    # write-out all pickles again since object attributes may have been updated
    pyorganism.write_pickle(genes, gene_file)
    pyorganism.write_pickle(products, product_file)
    pyorganism.write_pickle(t_factors, tf_file)
    if version >= (7, 2):
        pyorganism.write_pickle(sigma_factors, sigma_file)
    pyorganism.write_pickle(conformations, conf_file)
    pyorganism.write_pickle(promoters, prom_file)
    pyorganism.write_pickle(operons, op_file)
    pyorganism.write_pickle(t_units, tu_file)
    pyorganism.write_pickle(interactions, inter_file)

if __name__ == "__main__":
    if len(sys.argv) < 3 or len(sys.argv) > 4:
        LOGGER.error("{0} <distribution directory or archive>"\
                " <output directory> [version]".format(sys.argv[0]))
    else:
        main(sys.argv[1:])

