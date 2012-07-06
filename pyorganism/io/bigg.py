#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
===============================
BiGG Database Metabolite Parser
===============================

:Authors:
    Moritz Emanuel Beber
:Date:
    2012-07-05
:Copyright:
    Copyright(c) 2012 Jacobs University of Bremen. All rights reserved.
:File:
    bigg.py
"""


__all__ = ["BiGGMetaboliteParser"]


import logging
import csv

from ..singletonmixin import Singleton
from .. import miscellaneous as misc
from ..metabolism.elements import SBMLCompound
from .generic import open_file


LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(misc.NullHandler())

OPTIONS = misc.OptionsManager.get_instance()


class BiGGMetaboliteParser(Singleton):

    def __init__(self, **kw_args):
        super(BiGGMetaboliteParser, self).__init__(**kw_args)
        self.compounds = set()
        self.mapping = dict()

    def __call__(self, filename, sep="\t", comment="#", mode="rb",
            encoding="utf-8", **kw_args):
        compounds = set()
        mapping = dict()
        kw_args["mode"] = mode
        kw_args["encoding"] = encoding
        with  open_file(filename, **kw_args) as (file_h, ext):
            dialect = csv.Sniffer().sniff(file_h.read(2048), delimiters=sep)
            file_h.seek(0)
            reader = csv.DictReader(file_h, dialect=dialect)
            for row in reader:
                try:
                    charge = int(row["charge"])
                except ValueError:
                    charge = None
                cmpd = SBMLCompound(unique_id=row["abbreviation"],
                        name=row["name"], formula=row["formula"],
                        kegg_id=row["kegg_id"], cas_id=row["cas_id"],
                        charge=charge)
                compounds.add(cmpd)
                if not cmpd.kegg_id is None:
                    mapping[cmpd.kegg_id] = cmpd.unique_id
        self.compounds.update(compounds)
        self.mapping.update(mapping)
        return (compounds, mapping)

