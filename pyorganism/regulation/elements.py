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
        "NucleoidAssociatedProtein"]


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

    def __contains__(self, name):
        if name == self.unique_id:
            return True
        elif name == self.name:
            return True
        # need substring test for bnumber for entries with additional info
        elif name in self.bnumber:
            return True
        elif self.synonyms and any(name in syn for syn in self.synonyms if syn):
            return True
        else:
            return False

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
        elif self.synonyms and any(name in syn for syn in self.synonyms if syn):
            return True
        elif name == self.go:
            return True
        else:
            return False

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
        elif not self.synonyms is None and any(name in syn for syn in\
                self.synonyms if syn):
            return True
        elif name == self.go:
            return True
        else:
            return False

    def print_info(self, stream=sys.stdout):
        print >> stream, "ECK12:", self.unique_id
        print >> stream, "name:", self.name
        print >> stream, "synonyms:", self.synonyms
        print >> stream, self.go


class TranscriptionFactor(Regulator):

    def __init__(self, unique_id="", **kw_args):
        super(TranscriptionFactor, self).__init__(unique_id=unique_id, **kw_args)


class SigmaFactor(Regulator):

    def __init__(self, unique_id="", **kw_args):
        super(SigmaFactor, self).__init__(unique_id=unique_id, **kw_args)


class NucleoidAssociatedProtein(Regulator):

    def __init__(self, unique_id="", **kw_args):
        super(NucleoidAssociatedProtein, self).__init__(unique_id=unique_id,
                **kw_args)


