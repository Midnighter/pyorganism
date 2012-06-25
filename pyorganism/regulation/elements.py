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


__all__ = ["Gene", "TranscriptionalRegulator"]


import logging

from .. import miscellaneous as misc
from ..base import UniqueBase


LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(misc.NullHandler())


class Gene(UniqueBase):

    def __init__(self, name="", eck12="", bnumber="", synonyms=None,
            position=None, strand=None, product=None,
            *args, **kw_args):
        super(Gene, self).__init__(name=name, *args, **kw_args)
        self.eck12 = eck12
        self.bnumber = bnumber
        self.synonyms = misc.convert(synonyms, list, list())
        self.position = position
        self.strand = strand
        self.product = product

    def __contains__(self, element):
        if element == self.name:
            return True
        elif element == self.bnumber:
            return True
        elif element in self.synonyms:
            return True
        elif element == self.eck12:
            return True
        else:
            return False


class TranscriptionalRegulator(UniqueBase):

    def __init__(self, name="", synonyms=None, go=None, *args, **kw_args):
        super(TranscriptionalRegulator, self).__init__(name=name, *args,
                **kw_args)
        self.synonyms = misc.convert(synonyms, list, list())
        self.go = go

    def __contains__(self, element):
        if element == self.name:
            return True
        elif element in self.synonyms:
            return True
        elif element == self.go:
            return True
        else:
            return False

