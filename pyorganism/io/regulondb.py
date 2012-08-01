#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
==================================
RegulonDB Parsing and Web Services
==================================

:Authors:
    Moritz Emanuel Beber
:Date:
    2012-06-03
:Copyright:
    Copyright(c) 2012 Jacobs University of Bremen. All rights reserved.
:File:
    regulondb.py
"""


import logging
import re
#import xml.etree.ElementTree as etree

#from StringIO import StringIO
#from xml.etree.ElementTree import ElementTree
from bs4 import BeautifulSoup
from .. import miscellaneous as misc
from .generic import open_file, read_tabular
from ..regulation import elements as elem
from ..regulation.networks import GRN


LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(misc.NullHandler())

NAPs = ["ECK120011186", "ECK120011224", "ECK120011229", "ECK120011235",
        "ECK120011294", "ECK120011345", "ECK120011383"]


################################################################################
# RegulonDB dumps in XML format are not actually XML conform. There are problems
# with case-sensitivity of tags, for example, <i>...</I>, non-XML conform
# breaks, i.e., <br> instead of <br />, and HTML entities, like &sigma;. In
# conclusion: Abandon XML parsers and use BeautifulSoup to parse them as XHTML.
# BeautifulSoup is unfortunately much slower but using a stricter XML parser is
# insane.
################################################################################


#def sanitize_xml(xml_string):
#    """
#    Cleans an XML string and returns a file-like object.
#    """
#    # clean <BR> tags
#    break_tag = re.compile(r"<[bB]{1}[rR]{1}>")
#    (xml_string, nsubs) = break_tag.subn("<BR />", xml_string)
#    italic_start_tag = re.compile(r"<[iI]{1}>")
#    (xml_string, nsubs) = italic_start_tag.subn("<I>", xml_string)
#    italic_end_tag = re.compile(r"</[iI]{1}>")
#    (xml_string, nsubs) = italic_end_tag.subn("</I>", xml_string)
##    notes = re.compile(r"<(PRODUCT|GENE)_NOTE>(.*?)</\1_NOTE>")
##    (xml_string, nsubs) = notes.subn("<\g<1>_NOTE></\g<1>_NOTE>", xml_string)
##    xml_string = xml_string.replace("&", "&amp;")
#    return StringIO(xml_string)

def read_genes(filename, sep="\t", comment="#", encoding=None, mode="rb",
        **kw_args):
    """
    Retrieve the gene locations from a RegulonDB flat file and construct the
    GPN using the `proximity_threshold`.

    Parameters
    ----------
    filename: str
        The path to the file containing the interaction information.
    gene_identifier: str (optional)
        Identifiers for the genes used in the network. Can be any of:
            * 'name' for the gene name in that organism
            * 'blattner' for the Blattner number
            * 'regulondb' for the unique identifier assigned by RegulonDB
    sep: str (optional)
        The column separator; not likely to change.
    comment: str (optional)
        The sign denoting a comment line.
    encoding: str (optional)
        The file encoding.
    mode: str (optional)
        The mode used to open a file, should always be read binary.
    """

    def parse_flat_file(file_handle):
        raise NotImplementedError
        for row in read_tabular(file_handle, sep=sep, comment=comment):
            try:
                gene = elem.Gene(unique_id=row[0], name=row[1], position_start=row[2],
                        position_end=row[3], strand=row[4], gc_content=row[6])
            except StandardError as err:
                LOGGER.warn(err)
                LOGGER.warn(row)
            else:
                genes.append(gene)

    def parse_xml_file(file_handle):
        tree = etree.parse(file_handle)
#        tree = ElementTree()
#        tree.parse(file_handle)
        for row in tree.iter("ROW"):
            start = row.findtext("GENE_POSLEFT")
            if not start:
                start = None
            end = row.findtext("GENE_POSRIGHT")
            if not end:
                end = None
            gc = row.findtext("GC_CONTENT")
            if not gc:
                gc = None
            gene = elem.Gene(unique_id=row.findtext("GENE_ID"),
                name=row.findtext("GENE_NAME"),
                position_start=start,
                position_end=end,
                strand=row.findtext("GENE_STRAND"),
                gc_content=gc)
            genes.append(gene)

    def parse_xhtml(file_handle):
        soup = BeautifulSoup(file_handle, "lxml")
        for row in soup.rowset.find_all(name="row", recursive=False):
            start = misc.convert(row.gene_posleft.string, unicode)
            if not start:
                start = None
            end = misc.convert(row.gene_posright.string, unicode)
            if not end:
                end = None
            gc = misc.convert(row.gc_content.string, unicode)
            if not gc:
                gc = None
            gene = elem.Gene(unique_id=misc.convert(row.gene_id.string, unicode),
                name=misc.convert(row.gene_name.string, unicode),
                position_start=start,
                position_end=end,
                strand=misc.convert(row.gene_strand.string, unicode),
                gc_content=gc)
            genes.append(gene)

    genes = list()
    # read information from the file
    kw_args["mode"] = mode
    kw_args["encoding"] = encoding
    with open_file(filename, **kw_args) as (file_h, ext):
        if ext == ".xml":
            parse_xhtml(file_h)
#            clean = sanitize_xml(file_h.read())
#            parse_xml_file(clean)
        else:
            parse_flat_file(file_h)
    return genes

def update_gene_synonyms(filename, sep="\t", comment="#", encoding=None,
        mode="rb", **kw_args):

    def parse_flat_file(file_handle):
        raise NotImplementedError

    def parse_xml_file(file_handle):
        tree = ElementTree()
        tree.parse(file_handle)
        for row in tree.iter("ROW"):
            try:
                gene = elem.Gene.get(row.findtext("OBJECT_ID"))
            except KeyError:
                continue
            synonym = row.findtext("OBJECT_SYNONYM_NAME")
            if bpattern.match(synonym):
                gene.bnumber = synonym
            else:
                gene.synonyms.add(synonym)

    def parse_xhtml(file_handle):
        soup = BeautifulSoup(file_handle, "lxml")
        for row in soup.rowset.find_all(name="row", recursive=False):
            try:
                gene = elem.Gene.get(misc.convert(row.object_id.string, unicode))
            except KeyError:
                continue
            synonym = misc.convert(row.object_synonym_name.string, unicode)
            if synonym and bpattern.match(synonym):
                if gene.bnumber or "obsolete" in synonym:
                    gene.synonyms.add(synonym)
                else:
                    gene.bnumber = synonym
            else:
                gene.synonyms.add(synonym)

    bpattern = re.compile(r"b\d{4}")
    # read information from the file
    kw_args["mode"] = mode
    kw_args["encoding"] = encoding
    with open_file(filename, **kw_args) as (file_h, ext):
        if ext == ".xml":
            parse_xhtml(file_h)
#            clean = sanitize_xml(file_h.read())
#            parse_xml_file(clean)
        else:
            parse_flat_file(file_h)

def update_gene_external(filename, sep="\t", comment="#", encoding=None,
        mode="rb", **kw_args):

    def parse_flat_file(file_handle):
        raise NotImplementedError

    def parse_xml_file(file_handle):
        raise NotImplementedError

    def parse_xhtml(file_handle):
        soup = BeautifulSoup(file_handle, "lxml")
        for row in soup.rowset.find_all(name="row", recursive=False):
            try:
                gene = elem.Gene.get(misc.convert(row.object_id.string, unicode))
            except KeyError:
                continue
            synonym = misc.convert(row.ext_reference_id.string, unicode)
            if synonym and bpattern.match(synonym):
                if gene.bnumber or "obsolete" in synonym:
                    gene.synonyms.add(synonym)
                else:
                    gene.bnumber = synonym
            else:
                gene.synonyms.add(synonym)

    bpattern = re.compile(r"b\d{4}")
    # read information from the file
    kw_args["mode"] = mode
    kw_args["encoding"] = encoding
    with open_file(filename, **kw_args) as (file_h, ext):
        if ext == ".xml":
            parse_xhtml(file_h)
#            clean = sanitize_xml(file_h.read())
#            parse_xml_file(clean)
        else:
            parse_flat_file(file_h)

def read_products(filename, sep="\t", comment="#", encoding=None,
        mode="rb", **kw_args):

    def parse_flat_file(file_handle):
        raise NotImplementedError

    def parse_xml_file(file_handle):
        tree = ElementTree()
        tree.parse(file_handle)
        for row in tree.iter("ROW"):
            # typo in RegulonDB files for molecular weight
            mw = row.findtext("MOLECULAR_WEIGTH")
            if not mw:
                mw = None
            ip = row.findtext("ISOELECTRIC_POINT")
            if not ip:
                ip = None
            product = elem.Product(unique_id=row.findtext("PRODUCT_ID"),
                name=row.findtext("PRODUCT_NAME"),
                molecular_weight=mw,
                isoelectric_point=ip)
            products.append(product)

    def parse_xhtml(file_handle):
        soup = BeautifulSoup(file_handle, "lxml")
        for row in soup.rowset.find_all(name="row", recursive=False):
            # typo in RegulonDB files for molecular weight
            mw = misc.convert(row.molecular_weigth.string, unicode)
            if not mw:
                mw = None
            ip = misc.convert(row.isoelectric_point.string, unicode)
            if not ip:
                ip = None
            product = elem.Product(unique_id=misc.convert(row.product_id.string, unicode),
                name=misc.convert(row.product_name.string, unicode),
                molecular_weight=mw,
                isoelectric_point=ip)
            products.append(product)

    products = list()
    # read information from the file
    kw_args["mode"] = mode
    kw_args["encoding"] = encoding
    with open_file(filename, **kw_args) as (file_h, ext):
        if ext == ".xml":
            parse_xhtml(file_h)
#            clean = sanitize_xml(file_h.read())
#            parse_xml_file(clean)
        else:
            parse_flat_file(file_h)
    return products

def update_product_synonyms(filename, sep="\t", comment="#", encoding=None,
        mode="rb", **kw_args):

    def parse_flat_file(file_handle):
        raise NotImplementedError

    def parse_xml_file(file_handle):
        tree = ElementTree()
        tree.parse(file_handle)
        for row in tree.iter("ROW"):
            try:
                product = elem.Product.get(row.findtext("OBJECT_ID"))
            except KeyError:
                continue
            product.synonyms.add(row.findtext("OBJECT_SYNONYM_NAME"))

    def parse_xhtml(file_handle):
        soup = BeautifulSoup(file_handle, "lxml")
        for row in soup.rowset.find_all(name="row", recursive=False):
            try:
                product = elem.Product.get(misc.convert(row.object_id.string, unicode))
            except KeyError:
                continue
            product.synonyms.add(misc.convert(row.object_synonym_name.string, unicode))

    # read information from the file
    kw_args["mode"] = mode
    kw_args["encoding"] = encoding
    with open_file(filename, **kw_args) as (file_h, ext):
        if ext == ".xml":
            parse_xhtml(file_h)
#            clean = sanitize_xml(file_h.read())
#            parse_xml_file(clean)
        else:
            parse_flat_file(file_h)

def update_product_external(filename, sep="\t", comment="#", encoding=None,
        mode="rb", **kw_args):

    def parse_flat_file(file_handle):
        raise NotImplementedError

    def parse_xml_file(file_handle):
        raise NotImplementedError

    def parse_xhtml(file_handle):
        soup = BeautifulSoup(file_handle, "lxml")
        for row in soup.rowset.find_all(name="row", recursive=False):
            try:
                product = elem.Product.get(misc.convert(row.object_id.string, unicode))
            except KeyError:
                continue
            synonym = misc.convert(row.ext_reference_id.string, unicode)
            if synonym.startswith("GO:"):
                product.go = synonym
            else:
                product.synonyms.add(misc.convert(row.ext_reference_id.string, unicode))

    # read information from the file
    kw_args["mode"] = mode
    kw_args["encoding"] = encoding
    with open_file(filename, **kw_args) as (file_h, ext):
        if ext == ".xml":
            parse_xhtml(file_h)
#            clean = sanitize_xml(file_h.read())
#            parse_xml_file(clean)
        else:
            parse_flat_file(file_h)

def link_gene_product(filename, sep="\t", comment="#", encoding=None,
        mode="rb", **kw_args):

    def parse_flat_file(file_handle):
        raise NotImplementedError

    def parse_xml_file(file_handle):
        tree = ElementTree()
        tree.parse(file_handle)
        for row in tree.iter("ROW"):
            try:
                eck12 = row.findtext("GENE_ID")
                gene = elem.Gene.get(eck12)
            except KeyError:
                LOGGER.warn("unknown gene: {0}".format(eck12))
                continue
            try:
                eck12 = row.findtext("PRODUCT_ID")
                gene.product = elem.Product.get(eck12)
            except KeyError:
                LOGGER.warn("unknown product: {0}".format(eck12))
                continue
            if gene.product.coded_from:
                LOGGER.warn("product '{0}' already has coding gene '{1}'".format(
                        gene.product.unique_id,
                        gene.product.coded_from.unique_id))
            gene.product.coded_from = gene

    def parse_xhtml(file_handle):
        soup = BeautifulSoup(file_handle, "lxml")
        for row in soup.rowset.find_all(name="row", recursive=False):
            try:
                eck12 = misc.convert(row.gene_id.string, unicode)
                gene = elem.Gene.get(eck12)
            except KeyError:
                LOGGER.warn("unknown gene: {0}".format(eck12))
                continue
            try:
                eck12 = misc.convert(row.product_id.string, unicode)
                gene.product = elem.Product.get(eck12)
            except KeyError:
                LOGGER.warn("unknown product: {0}".format(eck12))
                continue
            if gene.product.coded_from:
                LOGGER.warn("product '{0}' already has coding gene '{1}'".format(
                        gene.product.unique_id,
                        gene.product.coded_from.unique_id))
            gene.product.coded_from = gene

    # read information from the file
    kw_args["mode"] = mode
    kw_args["encoding"] = encoding
    with open_file(filename, **kw_args) as (file_h, ext):
        if ext == ".xml":
            parse_xhtml(file_h)
#            clean = sanitize_xml(file_h.read())
#            parse_xml_file(clean)
        else:
            parse_flat_file(file_h)

def read_gene_regulation(filename, sep="\t", comment="#", encoding=None,
        mode="rb", **kw_args):

    def parse_flat_file(file_handle):
        raise NotImplementedError

    def parse_xml_file(file_handle):
        tree = ElementTree()
        tree.parse(file_handle)
        for row in tree.iter("ROW"):
            try:
                eck12 = row.findtext("GENE_ID_REGULATOR")
                gene_u = elem.Gene.get(eck12)
            except KeyError:
                LOGGER.warn("unknown gene: {0}".format(eck12))
                continue
            try:
                eck12 = row.findtext("GENE_ID_REGULATED")
                gene_v = elem.Gene.get(eck12)
            except KeyError:
                LOGGER.warn("unknown gene: {0}".format(eck12))
                continue
            interaction = row.findtext("GENEREGULATION_FUNCTION")
            if interaction == "repressor":
                trn.add_edge(gene_u, gene_v, interaction=-1)
            elif interaction == "activator":
                trn.add_edge(gene_u, gene_v, interaction=1)
            elif interaction == "unknown":
                trn.add_edge(gene_u, gene_v, interaction=0)
            elif interaction == "dual":
                trn.add_edge(gene_u, gene_v, interaction=2)
            else:
                LOGGER.warn("unknown regulatory interaction: {0}".format(interaction))

    def parse_xhtml(file_handle):
        soup = BeautifulSoup(file_handle, "lxml")
        for row in soup.rowset.find_all(name="row", recursive=False):
            try:
                eck12 = misc.convert(row.gene_id_regulator.string, unicode)
                gene_u = elem.Gene.get(eck12)
            except KeyError:
                LOGGER.warn("unknown gene: {0}".format(eck12))
                continue
            try:
                eck12 = misc.convert(row.gene_id_regulated.string, unicode)
                gene_v = elem.Gene.get(eck12)
            except KeyError:
                LOGGER.warn("unknown gene: {0}".format(eck12))
                continue
            interaction = misc.convert(row.generegulation_function.string, unicode)
            if interaction == "repressor":
                trn.add_edge(gene_u, gene_v, interaction=-1)
            elif interaction == "activator":
                trn.add_edge(gene_u, gene_v, interaction=1)
            elif interaction == "unknown":
                trn.add_edge(gene_u, gene_v, interaction=0)
            elif interaction == "dual":
                trn.add_edge(gene_u, gene_v, interaction=2)
            else:
                LOGGER.warn("unknown regulatory interaction: {0}".format(interaction))

    trn = GRN(name="E. coli Gene Regulatory Network")
    # read information from the file
    kw_args["mode"] = mode
    kw_args["encoding"] = encoding
    with open_file(filename, **kw_args) as (file_h, ext):
        if ext == ".xml":
            parse_xhtml(file_h)
#            clean = sanitize_xml(file_h.read())
#            parse_xml_file(clean)
        else:
            parse_flat_file(file_h)
    return trn

def read_sigma_factors(filename, sep="\t", comment="#", encoding=None,
        mode="rb", **kw_args):

    def parse_flat_file(file_handle):
        raise NotImplementedError

    def parse_xml_file(file_handle):
        tree = ElementTree()
        tree.parse(file_handle)
        for row in tree.iter("ROW"):
            sigma_factor = elem.SigmaFactor(unique_id=row.findtext("SIGMA_ID"),
                    name=row.findtext("SIGMA_NAME"))
            sigma_factor.synonyms.add(row.findtext("SIGMA_SYNONYMS"))
            try:
                eck12 = row.findtext("SIGMA_GENE_ID")
                gene = elem.Gene.get(eck12)
            except KeyError:
                LOGGER.warn("unknown gene: {0}".format(eck12))
                continue
            sigma_factor.made_from.add(gene.product)
            gene.regulatory_product = sigma_factor
            sigma_factor.coded_from.add(gene)
            s_factors.append(sigma_factor)

    def parse_xhtml(file_handle):
        soup = BeautifulSoup(file_handle, "lxml")
        for row in soup.rowset.find_all(name="row", recursive=False):
            sigma_factor = elem.SigmaFactor(unique_id=misc.convert(row.sigma_id.string, unicode),
                    name=misc.convert(row.sigma_name.string, unicode))
            sigma_factor.synonyms.add(misc.convert(row.sigma_synonyms.string, unicode))
            try:
                eck12 = misc.convert(row.sigma_gene_id.string, unicode)
                gene = elem.Gene.get(eck12)
            except KeyError:
                LOGGER.warn("unknown gene: {0}".format(eck12))
                continue
            sigma_factor.made_from.add(gene.product)
            gene.regulatory_product = sigma_factor
            sigma_factor.coded_from.add(gene)
            s_factors.append(sigma_factor)

    s_factors = list()
    # read information from the file
    kw_args["mode"] = mode
    kw_args["encoding"] = encoding
    with open_file(filename, **kw_args) as (file_h, ext):
        if ext == ".xml":
            parse_xhtml(file_h)
#            clean = sanitize_xml(file_h.read())
#            parse_xml_file(clean)
        else:
            parse_flat_file(file_h)
    return s_factors

def read_transcription_factors(filename, sep="\t", comment="#", encoding=None,
        mode="rb", **kw_args):

    def parse_flat_file(file_handle):
        raise NotImplementedError

    def parse_xml_file(file_handle):
        tree = ElementTree()
        tree.parse(file_handle)
        for row in tree.iter("ROW"):
            transcription_factor = elem.TranscriptionFactor(
                    unique_id=row.findtext("TRANSCRIPTION_FACTOR_ID"),
                    name=row.findtext("TRANSCRIPTION_FACTOR_NAME"))
            t_factors.append(transcription_factor)

    def parse_xhtml(file_handle):
        soup = BeautifulSoup(file_handle, "lxml")
        for row in soup.rowset.find_all(name="row", recursive=False):
            transcription_factor = elem.TranscriptionFactor(
                    unique_id=misc.convert(row.transcription_factor_id.string, unicode),
                    name=misc.convert(row.transcription_factor_name.string, unicode))
            t_factors.append(transcription_factor)

    t_factors = list()
    # read information from the file
    kw_args["mode"] = mode
    kw_args["encoding"] = encoding
    with open_file(filename, **kw_args) as (file_h, ext):
        if ext == ".xml":
            parse_xhtml(file_h)
#            clean = sanitize_xml(file_h.read())
#            parse_xml_file(clean)
        else:
            parse_flat_file(file_h)
    return t_factors

def update_transcription_factor_synonyms(filename, sep="\t", comment="#", encoding=None,
        mode="rb", **kw_args):

    def parse_flat_file(file_handle):
        raise NotImplementedError

    def parse_xml_file(file_handle):
        tree = ElementTree()
        tree.parse(file_handle)
        for row in tree.iter("ROW"):
            try:
                t_factor = elem.TranscriptionFactor.get(row.findtext("OBJECT_ID"))
            except KeyError:
                continue
            t_factor.synonyms.add(row.findtext("OBJECT_SYNONYM_NAME"))

    def parse_xhtml(file_handle):
        soup = BeautifulSoup(file_handle, "lxml")
        for row in soup.rowset.find_all(name="row", recursive=False):
            try:
                t_factor = elem.TranscriptionFactor.get(misc.convert(row.object_id.string, unicode))
            except KeyError:
                continue
            t_factor.synonyms.add(misc.convert(row.object_synonym_name.string, unicode))

    # read information from the file
    kw_args["mode"] = mode
    kw_args["encoding"] = encoding
    with open_file(filename, **kw_args) as (file_h, ext):
        if ext == ".xml":
            parse_xhtml(file_h)
#            clean = sanitize_xml(file_h.read())
#            parse_xml_file(clean)
        else:
            parse_flat_file(file_h)

def update_transcription_factor_external(filename, sep="\t", comment="#", encoding=None,
        mode="rb", **kw_args):

    def parse_flat_file(file_handle):
        raise NotImplementedError

    def parse_xml_file(file_handle):
        raise NotImplementedError

    def parse_xhtml(file_handle):
        soup = BeautifulSoup(file_handle, "lxml")
        for row in soup.rowset.find_all(name="row", recursive=False):
            try:
                t_factor = elem.TranscriptionFactor.get(misc.convert(row.object_id.string, unicode))
            except KeyError:
                continue
            t_factor.synonyms.add(misc.convert(row.ext_reference_id.string, unicode))

    # read information from the file
    kw_args["mode"] = mode
    kw_args["encoding"] = encoding
    with open_file(filename, **kw_args) as (file_h, ext):
        if ext == ".xml":
            parse_xhtml(file_h)
#            clean = sanitize_xml(file_h.read())
#            parse_xml_file(clean)
        else:
            parse_flat_file(file_h)

def link_gene_product_transcription_factor(filename, sep="\t", comment="#",
        encoding=None, mode="rb", **kw_args):

    def parse_flat_file(file_handle):
        raise NotImplementedError

    def parse_xml_file(file_handle):
        tree = ElementTree()
        tree.parse(file_handle)
        for row in tree.iter("ROW"):
            try:
                eck12 = row.findtext("TRANSCRIPTION_FACTOR_ID")
                t_factor = elem.TranscriptionFactor.get(eck12)
            except KeyError:
                LOGGER.warn("unknown transcription factor: {0}".format(eck12))
                continue
            try:
                eck12 = row.findtext("PRODUCT_ID")
                product = elem.Product.get(eck12)
            except KeyError:
                LOGGER.warn("unknown product: {0}".format(eck12))
                continue
            t_factor.made_from.add(product)
            t_factor.coded_from.add(product.coded_from)
            product.coded_from.regulatory_product = t_factor

    def parse_xhtml(file_handle):
        soup = BeautifulSoup(file_handle, "lxml")
        for row in soup.rowset.find_all(name="row", recursive=False):
            try:
                eck12 = misc.convert(row.transcription_factor_id.string, unicode)
                t_factor = elem.TranscriptionFactor.get(eck12)
            except KeyError:
                LOGGER.warn("unknown transcription factor: {0}".format(eck12))
                continue
            try:
                eck12 = misc.convert(row.product_id.string, unicode)
                product = elem.Product.get(eck12)
            except KeyError:
                LOGGER.warn("unknown product: {0}".format(eck12))
                continue
            t_factor.made_from.add(product)
            t_factor.coded_from.add(product.coded_from)
            product.coded_from.regulatory_product = t_factor

    # read information from the file
    kw_args["mode"] = mode
    kw_args["encoding"] = encoding
    with open_file(filename, **kw_args) as (file_h, ext):
        if ext == ".xml":
            parse_xhtml(file_h)
#            clean = sanitize_xml(file_h.read())
#            parse_xml_file(clean)
        else:
            parse_flat_file(file_h)

def read_sigma_gene_network(genes, filename, sep="\t", comment="#",
        encoding=None, mode="rb", **kw_args):

    def find_element(name, group):
        for element in group:
            if name in element:
                return element

    links = list()
    # read information from the file
    kw_args["mode"] = mode
    kw_args["encoding"] = encoding
    with open_file(filename, **kw_args) as (file_h, ext):
        if ext != ".txt":
            raise NotImplementedError
        for row in read_tabular(file_h, sep=sep, comment=comment):
            if row[2] == "-":
                links.append((row[0], row[1], -1, row[3]))
            elif row[2] == "+":
                links.append((row[0], row[1], 1, row[3]))
            elif row[2] == "?":
                links.append((row[0], row[1], 0, row[3]))
            elif row[2] == "+-":
                links.append((row[0], row[1], 2, row[3]))
            else:
                LOGGER.warn(row)
    # map sigma factor names to objects
    sigma_factors = [gene.regulatory_product for gene in genes\
            if type(gene.regulatory_product) == elem.SigmaFactor]
    sigma_names = set([quad[0] for quad in links])
    name2sigma = dict()
    for name in sigma_names:
        name2sigma[name] = find_element(name, sigma_factors)
    for (name, tar) in name2sigma.iteritems():
        if tar is None:
            LOGGER.warn("sigma factor '{0}' not found".format(name))
    # map gene names to objects
    gene_names = set([quad[1] for quad in links])
    name2gene = dict()
    for name in gene_names:
        name2gene[name] = find_element(name, genes)
    for (name, tar) in name2gene.iteritems():
        if tar is None:
            LOGGER.warn("gene '{0}' not found".format(name))
    # now apply the maps to the interactions
    interactions = [(name2sigma[quad[0]], name2gene[quad[1]], quad[2])\
            for quad in links]
    return interactions

#def parse_regulondb_entry(eck12, dscrptn):
#    dscrptn = dscrptn.split("\n")
#    gene = False
#    product = False
#    info = dict()
#    prod_info = dict()
#    blattner = re.compile(r"b\d+")
#    location = re.compile(r"\d+")
#    phantom = 0
#    for line in dscrptn:
#        if line.startswith("GENE"):
#            gene = True
#        elif line.startswith("PRODUCT"):
#            gene = False
#            product = True
#        elif gene and line.startswith("\t"):
#            line = line.strip()
#            tmp = line.split(":")
#            if tmp[0].startswith("Name"):
#                name = tmp[1].strip()
#                if name == "Phantom Gene":
#                    phantom += 1
#                    info["name"] = "%s %d" % (name, phantom)
#                else:
#                    info["name"] = name
#            elif tmp[0].startswith("Synonym"):
#                info["synonyms"] = list()
#                for name in tmp[1].split(","):
#                    name = name.strip()
#                    if blattner.match(name):
#                        info["bnumber"] = name
#                    else:
#                        info["synonyms"].append(name)
#            elif tmp[0].startswith("Map"):
#                info["position"] = location.findall(tmp[1])
#                info["position"] = map(int, info["position"])
#                info["position"] = tuple(info["position"])
#            elif tmp[0].startswith("Strand"):
#                info["strand"] = tmp[1].strip()
#        elif product and line.startswith("\t"):
#            line = line.strip()
#            tmp = line.split(":")
#            if tmp[0].startswith("Name"):
#                prod_info["name"] = tmp[1].strip()
#            elif tmp[0].startswith("Synonym"):
#                prod_info["synonyms"] = list()
#                for name in tmp[1].split(","):
#                    name = name.strip()
#                    prod_info["synonyms"].append(name)
#            elif tmp[0].startswith("GO"):
#                prod_info["go"] = int(tmp[2])
#        else:
#            gene = False
#            product = False
#            break
#    tf = TranscriptionFactor(**prod_info)
#    return Gene(eck12=eck12, product=tf, **info)

def lookup_gene_information(eck12_nums, num_threads=20,
        wsdl="http://regulondb.ccg.unam.mx/webservices/Gene.jws?wsdl"):
    """
    A threaded method that extracts gene information from RegulonDB.

    Parameters
    ----------
    eck12_nums: list
        A collection of ECK12 identifiers used by RegulonDB whose information
        should be retrieved.
    wsdl: str (optional)
        URL of the RegulonDB WSDL server.
    num_threads: int
        The number of desired simultaneous connections to the WSDL server.

    Notes
    -----
    Requires SOAPpy and an active internet connection.
    """

    from .wsdl import ThreadedWSDLFetcher
    from Queue import Queue
    # use a threaded approach to server querying
    tasks = Queue()
    for i in range(num_threads):
        thrd = ThreadedWSDLFetcher(tasks, wsdl)
        thrd.start()
    # get all gene descriptions from RegulonDB
    descriptions = list()
    for gene in eck12_nums:
        tasks.put(("getGene", gene, descriptions))
    tasks.join()
    return descriptions

