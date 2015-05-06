# -*- coding: utf-8 -*-


from __future__ import (absolute_import, unicode_literals)


"""
===============================
BiGG Database Compound Parser
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


__all__ = ["BiGGCompoundParser"]


import logging
import csv

from ..singletonmixin import Singleton
from .. import miscellaneous as misc
from ..metabolism.elements import SBMLCompound
from .generic import open_file


LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(misc.NullHandler())

OPTIONS = misc.OptionsManager.get_instance()


class BiGGCompoundParser(Singleton):
    """
    Parsing of compound information lists retrieved from the BiGG[1]_ database.

    Compounds parsed from a file are stored in the `compounds` attribute. This
    can be used to parse multiple files and retrieve all compounds at once.

    References
    ----------
    1.. Schellenberger, J., Park, J., Conrad, T., Palsson, B., 2010.
        BiGG: a Biochemical Genetic and Genomic knowledgebase of large scale
        metabolic reconstructions. BMC Bioinformatics 11, 213.
    """

    def __init__(self, **kw_args):
        super(BiGGCompoundParser, self).__init__(**kw_args)
        self.compounds = set()

    def __call__(self, filename, sep="\t", mode="rb", encoding="utf-8",
            **kw_args):
        compounds = set()
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
        self.compounds.update(compounds)
        return compounds

