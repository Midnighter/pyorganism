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
#from glob import glob

import pyorganism
import pyorganism.regulation as pyreg
#import pyorganism.metabolism as pymet
import pyorganism.io.regulondb as regdb


#class FileHandleFactory(object):
#
#    def __init__(self, base, **kw_args):
#        super(FileHandleFactory, self).__init__(**kw_args)
#        if os.path.isdir(base):
#            self.__call__ = self._filesystem
#        elif base.endswith(".zip"):
#            self.__call__ = self._zip_archive
#        elif base.endswith(".tar"):
#            self.__call__ = self._tar_archive
#
#    def __call__(self, filename):
#        pass
#
#    def _filesystem(self, filename):
#
#    def _zip_archive(self, filename):
#
#    def _tar_archive(self, filename):


def compile_genes(base_dir):
    LOGGER.info("Compiling genes...")
    genes = regdb.read_genes(os.path.join(base_dir, "gene.xml"))
    LOGGER.info("Found {0:d}".format(len(genes)))
    LOGGER.info("Compiling additional information...")
    regdb.update_gene_synonyms(os.path.join(base_dir, "object_synonym.xml"))
    regdb.update_gene_external(os.path.join(base_dir, "object_external_db_link.xml"))
    norm = float(len(genes))
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
    LOGGER.info("Found {0:d}".format(len(products)))
    LOGGER.info("Compiling additional information...")
    regdb.update_product_synonyms(os.path.join(base_dir, "object_synonym.xml"))
    regdb.update_product_external(os.path.join(base_dir, "object_external_db_link.xml"))
    norm = float(len(products))
    num = sum(1 for elem in products if elem.name)
    LOGGER.info("Found {0:d} products with names ({1:.2%})".format(num, num / norm))
    num = sum(1 for elem in products if elem.synonyms)
    LOGGER.info("Found {0:d} products with synonyms ({1:.2%})".format(num, num / norm))
    num = sum(1 for elem in products if elem.go)
    LOGGER.info("Found {0:d} products with GO ID ({1:.2%})".format(num, num / norm))
    return products

def compile_transcription_factors(base_dir):
    LOGGER.info("Compiling transcription factors...")
    t_factors = regdb.read_transcription_factors(os.path.join(base_dir, "transcription_factor.xml"))
    LOGGER.info("Found {0:d}".format(len(t_factors)))
    regdb.update_transcription_factor_synonyms(os.path.join(base_dir, "object_synonym.xml"))
    regdb.update_transcription_factor_external(os.path.join(base_dir, "object_external_db_link.xml"))
    norm = float(len(t_factors))
    num = sum(1 for elem in t_factors if elem.name)
    LOGGER.info("Found {0:d} transcription factors with names ({1:.2%})".format(num, num / norm))
    num = sum(1 for elem in t_factors if elem.synonyms)
    LOGGER.info("Found {0:d} transcription factors with synonyms ({1:.2%})".format(num, num / norm))
    return t_factors

def compile_sigma_factors(base_dir):
    LOGGER.info("Compiling sigma factors...")
    sigma_factors = regdb.read_sigma_factors(os.path.join(base_dir, "sigma_tmp.xml"))
    return sigma_factors

def compile_regulation(base_dir):
    LOGGER.info("Compiling regulatory interactions...")
    interactions = regdb.read_regulatory_interactions(os.path.join(base_dir, "regulatory_interaction.xml"))
    conformations = regdb.read_conformations(os.path.join(base_dir, "conformation.xml"))
    t_units = regdb.read_transcription_units(os.path.join(base_dir, "transcription_unit.xml"))
    tu_objects = regdb.read_tu_objects(os.path.join(base_dir, "tu_objects_tmp.xml"))

def main(argv):
    if not os.path.exists(argv[0]):
        sys.exit(1)
    if not os.path.exists(argv[1]):
        os.makedirs(argv[1])
    version = float(os.path.basename(argv[0]))
    # compile genes and gene products
    genes = compile_genes(argv[0])
    products = compile_products(argv[0])
    regdb.link_gene_product(os.path.join(argv[0], "gene_product_link.xml"))
    if version >= 7.2:
        sigma_f = compile_sigma_factors(argv[0])
    if not version == 6.0:
        t_f = compile_transcription_factors(argv[0])
        regdb.link_gene_product_transcription_factor(os.path.join(argv[0], "product_tf_link.xml"))
        regdb.product_tf_link_control(t_f, products)
        num = sum(1 for gene in genes if isinstance(gene.regulatory_product, pyreg.TranscriptionFactor))
        norm = float(len(genes))
        LOGGER.info("Found {0:d} genes that code for transcription factors({1:.2%})".format(num, num / norm))
    pyorganism.write_pickle(genes, os.path.join(argv[1], "genes.pkl"))

if __name__ == "__main__":
    if len(sys.argv) != 3:
        LOGGER.error("{0} <distribution directory or archive> <output directory>".format(sys.argv[0]))
    else:
        main(sys.argv[1:])

