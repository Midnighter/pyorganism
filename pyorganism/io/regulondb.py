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
import itertools

import numpy

from operator import attrgetter
from collections import defaultdict
from bs4 import BeautifulSoup
from .. import miscellaneous as misc
from ..errors import PyOrganismError
from .generic import open_file, read_tabular
from ..regulation import elements as elem
from ..regulation.networks import GRN, TRN


LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(misc.NullHandler())

BPATTERN = re.compile(r"b|B\d{4}")
PARSER = "lxml"
FUNCTIONS = {"repressor": -1, "activator": 1, "unknown": 0, "dual": 2,
        "-": -1, "+": 1, "?": 0, "+-": 2}
NAPs = frozenset(["ECK120011186", "ECK120011224", "ECK120011229", "ECK120011235",
        "ECK120011294", "ECK120011345", "ECK120011383"])


def find_element(name, group):
    for element in group:
        if name in element:
            return element

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

    def parse_xhtml(file_handle):
        soup = BeautifulSoup(file_handle, PARSER)
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
            name = misc.convert(row.gene_name.string, unicode)
            if name and bpattern.match(name):
                bnumber = name
            else:
                bnumber = None

            gene = elem.Gene(unique_id=misc.convert(row.gene_id.string, unicode),
                name=name,
                bnumber=bnumber,
                position_start=start,
                position_end=end,
                strand=misc.convert(row.gene_strand.string, unicode),
                gc_content=gc)
            genes.append(gene)

    genes = list()
    bpattern = BPATTERN
    # read information from the file
    kw_args["mode"] = mode
    kw_args["encoding"] = encoding
    with open_file(filename, **kw_args) as (file_h, ext):
        if ext == ".xml":
            parse_xhtml(file_h)
        else:
            parse_flat_file(file_h)
    return genes

def update_gene_synonyms(filename, sep="\t", comment="#", encoding=None,
        mode="rb", **kw_args):

    def parse_flat_file(file_handle):
        raise NotImplementedError

    def parse_xhtml(file_handle):
        soup = BeautifulSoup(file_handle, PARSER)
        for row in soup.rowset.find_all(name="row", recursive=False):
            obj_id = misc.convert(row.object_id.string, unicode)
            if obj_id not in elem.Gene:
                continue
            gene = elem.Gene[obj_id]
            synonym = misc.convert(row.object_synonym_name.string, unicode)
            if synonym is not None:
                if bpattern.match(synonym):
                    if gene.bnumber or "obsolete" in synonym:
                        gene.synonyms.add(synonym)
                    else:
                        gene.bnumber = synonym
                else:
                    gene.synonyms.add(synonym)

    bpattern = BPATTERN
    # read information from the file
    kw_args["mode"] = mode
    kw_args["encoding"] = encoding
    with open_file(filename, **kw_args) as (file_h, ext):
        if ext == ".xml":
            parse_xhtml(file_h)
        else:
            parse_flat_file(file_h)

def update_gene_external(filename, sep="\t", comment="#", encoding=None,
        mode="rb", **kw_args):

    def parse_flat_file(file_handle):
        raise NotImplementedError

    def parse_xhtml(file_handle):
        soup = BeautifulSoup(file_handle, PARSER)
        for row in soup.rowset.find_all(name="row", recursive=False):
            obj_id = misc.convert(row.object_id.string, unicode)
            if obj_id not in elem.Gene:
                continue
            gene = elem.Gene[obj_id]
            synonym = misc.convert(row.ext_reference_id.string, unicode)
            if synonym is not None:
                if bpattern.match(synonym):
                    if gene.bnumber or "obsolete" in synonym:
                        gene.synonyms.add(synonym)
                    else:
                        gene.bnumber = synonym
                else:
                    gene.synonyms.add(synonym)

    bpattern = BPATTERN
    # read information from the file
    kw_args["mode"] = mode
    kw_args["encoding"] = encoding
    with open_file(filename, **kw_args) as (file_h, ext):
        if ext == ".xml":
            parse_xhtml(file_h)
        else:
            parse_flat_file(file_h)

def read_products(filename, sep="\t", comment="#", encoding=None,
        mode="rb", **kw_args):

    def parse_flat_file(file_handle):
        raise NotImplementedError

    def parse_xhtml(file_handle):
        soup = BeautifulSoup(file_handle, PARSER)
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
        else:
            parse_flat_file(file_h)
    return products

def update_synonyms(filename, cls, sep="\t", comment="#", encoding=None,
        mode="rb", **kw_args):

    def parse_flat_file(file_handle):
        raise NotImplementedError

    def parse_xhtml(file_handle):
        soup = BeautifulSoup(file_handle, PARSER)
        for row in soup.rowset.find_all(name="row", recursive=False):
            obj_id = misc.convert(row.object_id.string, unicode)
            if obj_id not in cls:
                continue
            inst = cls[obj_id]
            synonym = misc.convert(row.object_synonym_name.string, unicode)
            if synonym is not None:
                inst.synonyms.add(synonym)

    # read information from the file
    kw_args["mode"] = mode
    kw_args["encoding"] = encoding
    with open_file(filename, **kw_args) as (file_h, ext):
        if ext == ".xml":
            parse_xhtml(file_h)
        else:
            parse_flat_file(file_h)

def update_external(filename, cls, sep="\t", comment="#", encoding=None,
        mode="rb", **kw_args):

    def parse_flat_file(file_handle):
        raise NotImplementedError

    def parse_xhtml(file_handle):
        soup = BeautifulSoup(file_handle, PARSER)
        for row in soup.rowset.find_all(name="row", recursive=False):
            obj_id = misc.convert(row.object_id.string, unicode)
            if obj_id not in cls:
                continue
            inst = cls[obj_id]
            synonym = misc.convert(row.ext_reference_id.string, unicode)
            if synonym is not None:
                inst.synonyms.add(synonym)

    # read information from the file
    kw_args["mode"] = mode
    kw_args["encoding"] = encoding
    with open_file(filename, **kw_args) as (file_h, ext):
        if ext == ".xml":
            parse_xhtml(file_h)
        else:
            parse_flat_file(file_h)

def link_gene_product(filename, sep="\t", comment="#", encoding=None,
        mode="rb", **kw_args):

    def parse_flat_file(file_handle):
        raise NotImplementedError

    def parse_xhtml(file_handle):
        soup = BeautifulSoup(file_handle, PARSER)
        for row in soup.rowset.find_all(name="row", recursive=False):
            try:
                eck12 = misc.convert(row.gene_id.string, unicode)
                gene = elem.Gene[eck12]
            except KeyError:
                LOGGER.warn("unknown gene: {0}".format(eck12))
                continue
            try:
                eck12 = misc.convert(row.product_id.string, unicode)
                gene.product = elem.Product[eck12]
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
        else:
            parse_flat_file(file_h)

def read_gene_regulation(filename, sep="\t", comment="#", encoding=None,
        mode="rb", **kw_args):

    def parse_flat_file(file_handle):
        raise NotImplementedError

    def parse_xhtml(file_handle):
        soup = BeautifulSoup(file_handle, PARSER)
        for row in soup.rowset.find_all(name="row", recursive=False):
            try:
                eck12 = misc.convert(row.gene_id_regulator.string, unicode)
                gene_u = elem.Gene[eck12]
            except KeyError:
                LOGGER.warn("unknown gene: {0}".format(eck12))
                continue
            try:
                eck12 = misc.convert(row.gene_id_regulated.string, unicode)
                gene_v = elem.Gene[eck12]
            except KeyError:
                LOGGER.warn("unknown gene: {0}".format(eck12))
                continue
            interaction = misc.convert(row.generegulation_function.string, unicode)
            grn.add_edge(gene_u, gene_v, key=FUNCTIONS[interaction])

    grn = GRN(name="Gene Regulatory Network")
    # read information from the file
    kw_args["mode"] = mode
    kw_args["encoding"] = encoding
    with open_file(filename, **kw_args) as (file_h, ext):
        if ext == ".xml":
            parse_xhtml(file_h)
        else:
            parse_flat_file(file_h)
    return grn

def read_sigma_factors(filename, sep="\t", comment="#", encoding=None,
        mode="rb", **kw_args):

    def parse_flat_file(file_handle):
        raise NotImplementedError

    def parse_xhtml(file_handle):
        soup = BeautifulSoup(file_handle, PARSER)
        for row in soup.rowset.find_all(name="row", recursive=False):
            sigma_factor = elem.SigmaFactor(unique_id=misc.convert(row.sigma_id.string, unicode),
                    name=misc.convert(row.sigma_name.string, unicode))
            sigma_factor.synonyms.add(misc.convert(row.sigma_synonyms.string, unicode))
            try:
                eck12 = misc.convert(row.sigma_gene_id.string, unicode)
                gene = elem.Gene[eck12]
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
        else:
            parse_flat_file(file_h)
    return s_factors

def read_transcription_factors(filename, sep="\t", comment="#", encoding=None,
        mode="rb", **kw_args):
    """
    TODO: transcription_factor.xml changes rapidly across versions, so depending
    on version we need to update various bits of information on the class
    instances.
    """

    def parse_flat_file(file_handle):
        raise NotImplementedError

    def parse_xhtml(file_handle):
        soup = BeautifulSoup(file_handle, PARSER)
        for row in soup.rowset.find_all(name="row", recursive=False):
            t_factor = elem.TranscriptionFactor(
                    unique_id=misc.convert(row.transcription_factor_id.string, unicode),
                    name=misc.convert(row.transcription_factor_name.string, unicode))
            t_factors.append(t_factor)

    # read information from the file
    kw_args["mode"] = mode
    kw_args["encoding"] = encoding
    t_factors = list()
    with open_file(filename, **kw_args) as (file_h, ext):
        if ext == ".xml":
            parse_xhtml(file_h)
        else:
            parse_flat_file(file_h)
    return t_factors

def update_product_transcription_factor_link(filename, sep="\t", comment="#",
        encoding=None, mode="rb", **kw_args):

    def parse_flat_file(file_handle):
        raise NotImplementedError

    def parse_xhtml(file_handle):
        soup = BeautifulSoup(file_handle, PARSER)
        for row in soup.rowset.find_all(name="row", recursive=False):
            prod_id = misc.convert(row.product_id.string, unicode)
            if prod_id not in elem.Product:
                LOGGER.warn("unknown product: {0}".format(prod_id))
                continue
            tf_id = misc.convert(row.transcription_factor_id.string, unicode)
            if prod_id not in elem.Product:
                LOGGER.warn("unknown transcription factor: {0}".format(prod_id))
                continue
            product = elem.Product[prod_id]
            t_factor = elem.TranscriptionFactor[tf_id]
            t_factor.made_from.add(product)
            t_factor.coded_from.add(product.coded_from)
            product.coded_from.regulatory_product = t_factor

    # read information from the file
    kw_args["mode"] = mode
    kw_args["encoding"] = encoding
    with open_file(filename, **kw_args) as (file_h, ext):
        if ext == ".xml":
            parse_xhtml(file_h)
        else:
            parse_flat_file(file_h)

def read_regulatory_interactions(filename, sep="\t", comment="#",
        encoding=None, mode="rb", **kw_args):
    """
    Extract regulatory interactions from relationships in RegulonDB files. Each
    interaction is given by a conformation ID and a promoter ID and the type of
    interaction.
    """

    def parse_flat_file(file_handle):
        raise NotImplementedError

    def parse_xhtml(file_handle):
        soup = BeautifulSoup(file_handle, PARSER)
        for row in soup.rowset.find_all(name="row", recursive=False):
            u = misc.convert(row.conformation_id.string, unicode)
            v = misc.convert(row.promoter_id.string, unicode)
            k = FUNCTIONS[misc.convert(row.ri_function.string, unicode)]
            interactions.append((u, v, k))

    # read information from the file
    kw_args["mode"] = mode
    kw_args["encoding"] = encoding
    interactions = list()
    with open_file(filename, **kw_args) as (file_h, ext):
        if ext == ".xml":
            parse_xhtml(file_h)
        else:
            parse_flat_file(file_h)
    return interactions

def read_transcription_units(filename, sep="\t", comment="#", encoding=None,
        mode="rb", **kw_args):
    """
    Extract pairs of transcription unit and promoter IDs.

    Returns
    -------
    A map between promoter IDs and the various transcription units they occur
    in.
    """

    def parse_flat_file(file_handle):
        raise NotImplementedError

    def parse_xhtml(file_handle):
        soup = BeautifulSoup(file_handle, PARSER)
        for row in soup.rowset.find_all(name="row", recursive=False):
            tu_id = misc.convert(row.transcription_unit_id.string, unicode)
            prom = misc.convert(row.promoter_id.string, unicode)
            op = misc.convert(row.operon_id.string, unicode)
            try:
                if prom is not None:
                    prom = elem.Promoter[prom]
                if op is not None:
                    op = elem.Operon[op]
            except KeyError:
                LOGGER.error("unknown promoter %s or operon %s", prom, op)
                LOGGER.error("Please parse operon and promoter information before\
                        parsing transcription units.")
                continue
            tu = elem.TranscriptionUnit(unique_id=tu_id,
                    name=misc.convert(row.transcription_unit_id.string, unicode),
                    promoter=prom, operon=op)
            units.append(tu)

    # read information from the file
    kw_args["mode"] = mode
    kw_args["encoding"] = encoding
    units = list()
    with open_file(filename, **kw_args) as (file_h, ext):
        if ext == ".xml":
            parse_xhtml(file_h)
        else:
            parse_flat_file(file_h)
    return units

def link_tu_genes(filename, sep="\t", comment="#", encoding=None,
        mode="rb", **kw_args):
    """
    Link transcription units (TU) with their corresponding genes.

    Parameters
    ----------
    filename: str
        Relative or absolute path to file that contains the  RegulonDB information.

    Returns
    -------
    A dict with keys that are TU IDs and values that are gene IDs contained
    within that TU.
    """

    def parse_flat_file(file_handle):
        raise NotImplementedError

    def parse_xhtml(file_handle):
        soup = BeautifulSoup(file_handle, PARSER)
        for row in soup.rowset.find_all(name="row", recursive=False):
            tu_id = misc.convert(row.transcription_unit_id.string, unicode)
            try:
                t_unit = elem.TranscriptionUnit[tu_id]
            except KeyError:
                LOGGER.error("unknown transcription unit '%s', please parse"\
                        " those first.", tu_id)
            gene_id = misc.convert(row.gene_id.string, unicode)
            try:
                gene = elem.Gene[gene_id]
            except KeyError:
                LOGGER.error("unknown gene '%s', please parse those first.", gene_id)
            t_unit.genes.append(gene)
            t_units.append(t_unit)

    # read information from the file
    kw_args["mode"] = mode
    kw_args["encoding"] = encoding
    t_units = list()
    with open_file(filename, **kw_args) as (file_h, ext):
        if ext == ".xml":
            parse_xhtml(file_h)
        else:
            parse_flat_file(file_h)
    begin = attrgetter("position")
    for tu in t_units:
        if all(gene.strand == "reverse" for gene in tu.genes):
            tu.genes.sort(key=begin, reverse=True)
        elif all(gene.strand == "forward" for gene in tu.genes):
            tu.genes.sort(key=begin)
        else:
            LOGGER.error("conflicting strand information in transcription unit"\
                    " '%s'", tu.unique_id)

def read_tu_objects(filename, sep="\t", comment="#", encoding=None,
        mode="rb", **kw_args):
    """
    Extract transcription unit (TU) information.

    Parameters
    ----------
    filename: str
        Relative or absolute path to file that contains the  RegulonDB information.

    Returns
    -------
    A dict with keys that are TU IDs and values that are gene IDs contained
    within that TU.
    """

    def parse_flat_file(file_handle):
        raise NotImplementedError

    def parse_xhtml(file_handle):
        soup = BeautifulSoup(file_handle, PARSER)
        for row in soup.rowset.find_all(name="row", recursive=False):
            if misc.convert(row.tu_object_class.string, unicode) == "GN":
                tu = misc.convert(row.transcription_unit_id.string, unicode)
                gene = misc.convert(row.tu_object_id.string, unicode)
                units[tu].append(gene)

    # read information from the file
    kw_args["mode"] = mode
    kw_args["encoding"] = encoding
    units = defaultdict(list)
    with open_file(filename, **kw_args) as (file_h, ext):
        if ext == ".xml":
            parse_xhtml(file_h)
        else:
            parse_flat_file(file_h)
    return dict(units)

def read_promoters(filename, sep="\t", comment="#", encoding=None,
        mode="rb", **kw_args):
    """
    Extract promoter information.

    Parameters
    ----------
    filename: str
        Relative or absolute path to file that contains the  RegulonDB information.

    Returns
    -------
    """

    def parse_flat_file(file_handle):
        raise NotImplementedError

    def parse_xhtml(file_handle):
        soup = BeautifulSoup(file_handle, PARSER)
        for row in soup.rowset.find_all(name="row", recursive=False):
            strand = misc.convert(row.promoter_strand.string, unicode)
            if not strand:
                strand = None
            prom = elem.Promoter(unique_id=misc.convert(row.promoter_id.string, unicode),
                    name=misc.convert(row.promoter_name.string, unicode),
                    strand=strand,
                    pos_1=misc.convert(row.pos_1.string, int),
                    sequence=misc.convert(row.promoter_sequence.string, unicode),
                    sigma_factor=misc.convert(row.sigma_factor.string,
                            unicode, unicode()).split(","),
                    note=misc.convert(row.promoter_note.string, unicode))
            promoters.append(prom)

    # read information from the file
    kw_args["mode"] = mode
    kw_args["encoding"] = encoding
    promoters = list()
    with open_file(filename, **kw_args) as (file_h, ext):
        if ext == ".xml":
            parse_xhtml(file_h)
        else:
            parse_flat_file(file_h)
    return promoters

def read_operons(filename, sep="\t", comment="#", encoding=None,
        mode="rb", **kw_args):
    """
    Extract operon information.

    Parameters
    ----------
    filename: str
        Relative or absolute path to file that contains the  RegulonDB information.

    Returns
    -------
    """

    def parse_flat_file(file_handle):
        raise NotImplementedError

    def parse_xhtml(file_handle):
        soup = BeautifulSoup(file_handle, PARSER)
        for row in soup.rowset.find_all(name="row", recursive=False):
            gene_left = misc.convert(row.firstgeneposleft.string, int)
            if not gene_left:
                gene_left = None
            gene_right = misc.convert(row.lastgeneposright.string, int)
            if not gene_right:
                gene_right = None
            reg_left = misc.convert(row.regulationposleft.string, int)
            if not reg_left:
                reg_left = None
            reg_right = misc.convert(row.regulationposright.string, int)
            if not reg_right:
                reg_right = None
            strand = misc.convert(row.operon_strand.string, unicode)
            if not strand:
                strand = None
            op = elem.Operon(unique_id=misc.convert(row.operon_id.string, unicode),
                    name=misc.convert(row.operon_name.string, unicode),
                    strand=strand,
                    gene_position_start=gene_left,
                    gene_position_end=gene_right,
                    regulation_position_start=reg_left,
                    regulation_position_end=reg_right)
            operons.append(op)

    # read information from the file
    kw_args["mode"] = mode
    kw_args["encoding"] = encoding
    operons = list()
    with open_file(filename, **kw_args) as (file_h, ext):
        if ext == ".xml":
            parse_xhtml(file_h)
        else:
            parse_flat_file(file_h)
    return operons

def update_operons(operons, promoters, genes, **kw_args):
    """
    Extract operon information.

    Parameters
    ----------

    Returns
    -------
    """
    promoters = list(prom for prom in promoters if prom.pos_1 is not None)
    prom_i = numpy.arange(len(promoters), dtype=int)
    prom_pos = numpy.array([prom.pos_1 for prom in promoters])
    genes = list(gene for gene in genes if gene.position is not None)
    gene_i = numpy.arange(len(genes), dtype=int)
    gene_start = numpy.array([gene.position_start for gene in genes])
    gene_end = numpy.array([gene.position_end for gene in genes])
    begin = attrgetter("position")
    for op in operons:
        mask = numpy.logical_and((prom_pos >= op.regulation_position_start),
                (prom_pos <= op.regulation_position_end))
        for prom in (promoters[i] for i in prom_i[mask]):
            if prom.strand == op.strand:
                op.promoters.add(prom)
        mask = numpy.logical_and((gene_start >= op.gene_position_start),
                (gene_end <= op.gene_position_end))
        for gene in (genes[i] for i in gene_i[mask]):
            if gene.strand == op.strand:
                op.genes.append(gene)
                gene.operons.add(op)
        if op.strand == "reverse":
            op.genes.sort(key=begin, reverse=True)
        else:
            op.genes.sort(key=begin)

def read_conformations(filename, sep="\t", comment="#", encoding=None,
        mode="rb", **kw_args):
    """
    Extract conformation information.

    Parameters
    ----------
    filename: str
        Relative or absolute path to file that contains the  RegulonDB information.

    Returns
    -------
    """

    def parse_flat_file(file_handle):
        raise NotImplementedError

    def parse_xhtml(file_handle):
        soup = BeautifulSoup(file_handle, PARSER)
        for row in soup.rowset.find_all(name="row", recursive=False):
            tf_id = misc.convert(row.transcription_factor_id.string, unicode)
            try:
                t_factor = elem.TranscriptionFactor[tf_id]
            except KeyError:
                LOGGER.error("unknown transcription factor %s", tf_id)
                LOGGER.error("Please parse transcription factor information before"\
                        " parsing conformations.")
                continue
            conf = elem.Conformation(
                    unique_id=misc.convert(row.conformation_id.string, unicode),
                    tf=t_factor,
                    state=misc.convert(row.final_state.string, unicode),
                    interaction=misc.convert(row.interaction_type.string, unicode),
                    conformation_type=misc.convert(row.conformation_type.string, unicode),
                    apo_holo=misc.convert(row.apo_holo_conformation.string, unicode)
            )
            t_factor.conformations.add(conf)
            conformations.append(conf)

    # read information from the file
    kw_args["mode"] = mode
    kw_args["encoding"] = encoding
    conformations = list()
    with open_file(filename, **kw_args) as (file_h, ext):
        if ext == ".xml":
            parse_xhtml(file_h)
        else:
            parse_flat_file(file_h)
    return conformations

def read_tf_gene_network(filename, name2tf, name2gene, sep="\t",
        encoding=None, mode="rb", comment="#", **kw_args):
    """
    Retrieve the TRN from a RegulonDB flat file.

    Parameters
    ----------
    filename: str
        The path to the file containing the interaction information.
    name2tf: dict
        The transcription factors in the parsed interactions are given by
        their name. A map between names and objects is required.
    name2gene: dict
        The genes in the parsed interactions are given by their name. A map
        between names and objects is required.
    sep: str (optional)
        The column separator; not likely to change.
    encoding: str (optional)
        The file encoding.
    mode: str (optional)
        The mode used to open a file, should always be read binary.
    """
    trn = TRN(name=filename)
    kw_args["encoding"] = encoding
    kw_args["mode"] = mode
    with open_file(filename, **kw_args) as (file_h, ext):
        for row in read_tabular(file_h, sep=sep, comment=comment):
            first = name2tf[row[0]]
            second = name2gene[row[1]]
            trn.add_edge(first, second, key=FUNCTIONS[row[2]], evidence=row[3])
    return trn

#def read_sigma_gene_network(genes, filename, sep="\t", comment="#",
#        encoding=None, mode="rb", **kw_args):
#
#    links = list()
#    # read information from the file
#    kw_args["mode"] = mode
#    kw_args["encoding"] = encoding
#    with open_file(filename, **kw_args) as (file_h, ext):
#        if ext != ".txt":
#            raise NotImplementedError
#        for row in read_tabular(file_h, sep=sep, comment=comment):
#            if row[2] == "-":
#                links.append((row[0], row[1], -1, row[3]))
#            elif row[2] == "+":
#                links.append((row[0], row[1], 1, row[3]))
#            elif row[2] == "?":
#                links.append((row[0], row[1], 0, row[3]))
#            elif row[2] == "+-":
#                links.append((row[0], row[1], 2, row[3]))
#            else:
#                LOGGER.warn(row)
#    # map sigma factor names to objects
#    sigma_factors = [gene.regulatory_product for gene in genes\
#            if type(gene.regulatory_product) == elem.SigmaFactor]
#    sigma_names = set([quad[0] for quad in links])
#    name2sigma = dict()
#    for name in sigma_names:
#        name2sigma[name] = find_element(name, sigma_factors)
#    for (name, tar) in name2sigma.iteritems():
#        if tar is None:
#            LOGGER.warn("sigma factor '{0}' not found".format(name))
#    # map gene names to objects
#    gene_names = set([quad[1] for quad in links])
#    name2gene = dict()
#    for name in gene_names:
#        name2gene[name] = find_element(name, genes)
#    for (name, tar) in name2gene.iteritems():
#        if tar is None:
#            LOGGER.warn("gene '{0}' not found".format(name))
#    # now apply the maps to the interactions
#    interactions = [(name2sigma[quad[0]], name2gene[quad[1]], quad[2])\
#            for quad in links]
#    return interactions

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

def draw_relationships(file_contents, emph=list(), ignore=["key_id_org"],
        title="", font_size=14.0, width=16.54, height=11.69):
    """
    Draw a graph with relationships between tags in RegulonDB XML files.

    Parameters
    ----------
    file_contents: dict
        A dictionary where keys are filenames and values are lists of tag names
        of the corresponding file.
    emph: list
        A list of tags that should be emphasized by colour.
    ignore: list
        A list of tags to ignore when drawing relationships.
    title: str
        A title description for the printed graph.

    Notes
    -----
    Requires pygraphviz.
    """
    if len(emph) > len(misc.BREWER_SET1):
        raise PyOrganismError("number of objects to be emphasized ({0:d}) is"\
                " greater than the number of colours available ({1:d})",
                len(emph), len(misc.BREWER_SET1))
    pgv = misc.load_module("pygraphviz")
    colour_choice = dict(itertools.izip(emph, misc.BREWER_SET1))
    graph = pgv.AGraph(name="RegulonDB File-Relationships", strict=True,
            directed=False, rankdir="TB")
    graph.graph_attr["labelloc"] = "t"
    graph.graph_attr["label"] = title
    graph.graph_attr["fontsize"] = font_size * 1.5
    graph.graph_attr["ranksep"] = "0.1 equally"
    graph.graph_attr["size"] = (width, height)
    graph.graph_attr["ratio"] = "compress"
    graph.node_attr["shape"] = "none"
    graph.node_attr["fontsize"] = font_size
    for (name, attrs) in file_contents.iteritems():
        label = ["<<TABLE BORDER=\"0\" CELLBORDER=\"1\" CELLSPACING=\"0\""\
                " CELLPADDING=\"4\">"]
        label.append("<TR><TD BGCOLOR=\"#A4A4A4\"><B>{0}</B></TD></TR>".format(name))
        for (i, attr) in enumerate(attrs):
            if attr in emph:
                label.append("<TR><TD PORT=\"f{0:d}\" BGCOLOR=\"{1}\">{2}</TD></TR>".format(i,
                            colour_choice[attr], attr))
            else:
                label.append("<TR><TD PORT=\"f{0:d}\">{1}</TD></TR>".format(i,
                    attr))
        label.append("</TABLE>>")
        graph.add_node(name, label="\n".join(label))
    nodes = file_contents.keys()
    for i in range(len(nodes) - 1):
        node_u = nodes[i]
        attr_u = file_contents[node_u]
        for j in range(i + 1, len(nodes)):
            node_v = nodes[j]
            attr_v = file_contents[node_v]
            shared = set(attr_u).intersection(set(attr_v))
            for attr in shared:
                if attr in ignore:
                    continue
                u = attr_u.index(attr)
                v = attr_v.index(attr)
                if attr in emph:
                    graph.add_edge(node_u, node_v,
                            tailport="f{0:d}".format(u), headport="f{0:d}".format(v),
                            color=colour_choice[attr])
                else:
                    graph.add_edge(node_u, node_v,
                            tailport="f{0:d}".format(u), headport="f{0:d}".format(v))
    sub_attr = dict()
    nodes = graph.nodes()
    nodes.sort(key=lambda n: graph.degree(n))
    maxi = nodes[-1: -len(nodes) / 4]
    nodes = nodes[:-len(nodes) / 4]
    zeros = [node for (node, deg) in graph.degree_iter() if deg == 0]
    for n in zeros:
        nodes.remove(n)
    graph.add_subgraph(maxi, name="input", rank="source", **sub_attr)
    graph.add_subgraph(nodes, name="middle", **sub_attr)
    graph.add_subgraph(zeros, rank="sink", **sub_attr)
    return graph

