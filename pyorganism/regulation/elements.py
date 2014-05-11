#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
==============================
PyOrganism Regulatory Elements
==============================

:Authors:
    Moritz Emanuel Beber
:Date:
    2012-06-08
:Copyright:
    Copyright(c) 2012 Jacobs University of Bremen. All rights reserved.
:File:
    elements.py
"""


__all__ = ["Gene", "Product", "TranscriptionFactor", "SigmaFactor",
        "NucleoidAssociatedProtein", "Promoter", "TranscriptionUnit", "Operon",
        "Conformation", "clear_memory"]

import sys
import logging

from .. import miscellaneous as misc
#from ..errors import PyOrganismError
from ..base import UniqueBase


LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(misc.NullHandler())


class Gene(UniqueBase):

    def __init__(self, unique_id="", name="", bnumber="", synonyms=None,
            position_start=None, position_end=None, strand=None, sequence=None,
            gc_content=None, product=None, regulatory_product=None, **kw_args):
        super(Gene, self).__init__(unique_id=unique_id, **kw_args)
        self.name = name
        self.bnumber = bnumber
        self.synonyms = misc.convert(synonyms, set, set())
        self.position_start = misc.convert(position_start, int)
        self.position_end = misc.convert(position_end, int)
        self.position = (self.position_start, self.position_end)
        self.strand = strand
        self.sequence = sequence
        self.gc_content = misc.convert(gc_content, float)
        self.product = product
        self.regulatory_product = regulatory_product
        self.operons = set()

    def __contains__(self, name):
        if name == self.unique_id:
            return True
        elif name == self.name:
            return True
        # need substring test for bnumber for entries with additional info
        elif name == self.bnumber:
            return True
        elif self.synonyms and any(name == syn for syn in self.synonyms if syn):
            return True
        else:
            return False

    def get_operons(self):
        return self.operons

    def print_info(self, stream=sys.stdout):
        print >> stream, "ECK12:", self.unique_id
        print >> stream, "name:", self.name
        print >> stream, "bnumber:", self.bnumber
        print >> stream, "synonyms:", self.synonyms
        print >> stream, "position:", self.position


class Product(UniqueBase):

    def __init__(self, unique_id="", name="", molecular_weight=None,
            isoelectric_point=None, synonyms=None, go=None, coded_from=None,
            **kw_args):
        super(Product, self).__init__(unique_id=unique_id, **kw_args)
        self.name = name
        self.molecular_weight = misc.convert(molecular_weight, float)
        self.isoelectric_point = misc.convert(isoelectric_point, float)
        self.synonyms = misc.convert(synonyms, set, set())
        self.go = go
        self.coded_from = coded_from

    def __contains__(self, name):
        if name == self.unique_id:
            return True
        elif name == self.name:
            return True
        elif self.synonyms and any(name == syn for syn in self.synonyms if syn):
            return True
        elif name == self.go:
            return True
        else:
            return False

    def get_operons(self):
        return set(op for gene in self.coded_from for op in gene.operons)

    def print_info(self, stream=sys.stdout):
        print >> stream, "ECK12:", self.unique_id
        print >> stream, "name:", self.name
        print >> stream, "synonyms:", self.synonyms
        print >> stream, self.go


class Regulator(UniqueBase):

    def __init__(self, unique_id="", name="", synonyms=None, go=None,
            coded_from=None, made_from=None, **kw_args):
        super(Regulator, self).__init__(unique_id=unique_id, **kw_args)
        self.name = name
        self.synonyms = misc.convert(synonyms, set, set())
        self.go = go
        self.coded_from = misc.convert(coded_from, set, set())
        self.made_from = misc.convert(made_from, set, set())

    def __contains__(self, name):
        if name == self.unique_id:
            return True
        elif name == self.name:
            return True
        elif not self.synonyms is None and any(name == syn for syn in\
                self.synonyms if syn):
            return True
        elif name == self.go:
            return True
        else:
            return False

    def get_operons(self):
        return set(op for gene in self.coded_from for op in gene.operons)

    def print_info(self, stream=sys.stdout):
        print >> stream, "ECK12:", self.unique_id
        print >> stream, "name:", self.name
        print >> stream, "synonyms:", self.synonyms
        print >> stream, self.go


class TranscriptionFactor(Regulator):

    def __init__(self, unique_id="", conformations=None, **kw_args):
        super(TranscriptionFactor, self).__init__(unique_id=unique_id, **kw_args)
        self.conformations = misc.convert(conformations, set, set())


class SigmaFactor(Regulator):

    def __init__(self, unique_id="", **kw_args):
        super(SigmaFactor, self).__init__(unique_id=unique_id, **kw_args)


class NucleoidAssociatedProtein(Regulator):

    def __init__(self, unique_id="", **kw_args):
        super(NucleoidAssociatedProtein, self).__init__(unique_id=unique_id,
                **kw_args)


class Conformation(UniqueBase):

    def __init__(self, unique_id="", name="", tf=None, state=None,
            conformation_type=None, interaction=None, apo_holo=None, **kw_args):
        super(Conformation, self).__init__(unique_id=unique_id, **kw_args)
        self.name = name
        self.t_factor = tf
        self.final_state = state
        self.type = conformation_type
        self.interaction = interaction
        self.apo_holo = apo_holo


class Promoter(UniqueBase):

    def __init__(self, unique_id="", name="", strand=None, pos_1=None,
            sequence=None, sigma_factor=None, note=None, **kw_args):
        super(Promoter, self).__init__(unique_id=unique_id,
                **kw_args)
        self.name = name
        self.strand = strand
        self.pos_1 = misc.convert(pos_1, int)
        self.sigma_factor = misc.convert(sigma_factor, list, list())
        self.sequence = sequence
        self.note = note

    def print_info(self, stream=sys.stdout):
        print >> stream, "ECK12:", self.unique_id
        print >> stream, "name:", self.name


class TranscriptionUnit(UniqueBase):

    def __init__(self, unique_id="", name="", promoter=None, operon=None,
            genes=None, **kw_args):
        super(TranscriptionUnit, self).__init__(unique_id=unique_id,
                **kw_args)
        self.name = name
        self.promoter = promoter
        self.operon = operon
        self.genes = misc.convert(genes, list, list())

    def __len__(self):
        return len(self.genes)

    def print_info(self, stream=sys.stdout):
        print >> stream, "ECK12:", self.unique_id
        print >> stream, "name:", self.name
        print >> stream, "Genes:", ", ".join([gene.name if gene.name else "?" for gene in self.genes])


class Operon(UniqueBase):

    def __init__(self, unique_id="", name="", strand=None, promoters=None, genes=None,
            gene_position_start=None, gene_position_end=None,
            regulation_position_start=None, regulation_position_end=None, **kw_args):
        super(Operon, self).__init__(unique_id=unique_id,
                **kw_args)
        self.name = name
        self.strand = strand
        self.gene_position_start = misc.convert(gene_position_start, int)
        self.gene_position_end = misc.convert(gene_position_end, int)
        self.regulation_position_start = misc.convert(regulation_position_start, int)
        self.regulation_position_end = misc.convert(regulation_position_end, int)
        self.promoters = misc.convert(promoters, set, set())
        self.genes = misc.convert(genes, list, list())

    def __len__(self):
        return len(self.genes)

    def print_info(self, stream=sys.stdout):
        print >> stream, "ECK12:", self.unique_id
        print >> stream, "name:", self.name
        print >> stream, "Genes:", ", ".join([gene.name if gene.name else "?" for gene in self.genes])


def clear_memory():
    Gene.clear()
    Product.clear()
    Regulator.clear()
    TranscriptionFactor.clear()
    SigmaFactor.clear()
    NucleoidAssociatedProtein.clear()
    Conformation.clear()
    Promoter.clear()
    TranscriptionUnit.clear()
    Operon.clear()

