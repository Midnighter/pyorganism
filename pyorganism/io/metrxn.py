# -*- coding: utf-8 -*-


from __future__ import (absolute_import, unicode_literals)


"""
=======================
MetRxn Database Parsers
=======================

:Authors:
    Moritz Emanuel Beber
:Date:
    2012-07-11
:Copyright:
    Copyright(c) 2012 Jacobs University of Bremen. All rights reserved.
:File:
    metrxn.py
"""


__all__ = ["MetRxnCompoundIDParser"]


import logging
import re

from ..singletonmixin import Singleton
from .. import miscellaneous as misc
from .generic import open_file, read_tabular


LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(misc.NullHandler())


class MetRxnCompoundIDParser(Singleton):
    """
    Parsing of compound identifier mappings retrieved from the MetRxn[1]_ database.

    References
    ----------
    1.. Kumar, A., Suthers, P., Maranas, C., 2012.
        "MetRxn: a knowledgebase of metabolites and reactions spanning metabolic
        models and databases". BMC Bioinformatics 13, 6.
    """

    def __init__(self, **kw_args):
        super(MetRxnCompoundIDParser, self).__init__(**kw_args)
        self.kegg_pattern = re.compile(r"C\d+")

    def __call__(self, filename, sep="\t", comment="#", mode="rb",
            encoding="utf-8", **kw_args):
        raise NotImplementedError
        mappings = dict()
        databases = None
        reverse = None
        kw_args["mode"] = mode
        kw_args["encoding"] = encoding
        with  open_file(filename, **kw_args) as (file_h, ext):
            for row in read_tabular(file_h, sep=sep, comment=comment):
                if not databases:
                    databases = tuple(db.strip() for db in\
                            row["Source names"].split(","))
                break
            reverse = (databases[1], databases[0])
            mappings[databases] = dict()
            mappings[reverse] = dict()
            for row in reader:
                ids = list()
                kegg_ids = list()
                # typo on database side
                tmp = [x.strip() for x in row["Metatbolite Id's"].split(",")]
                for cmpd in tmp:
                    if self.kegg_pattern.match(cmpd):
                        kegg_ids.append(cmpd)
                    else:
                        cmpd = cmpd.split("_")
                        ids.append(cmpd[1])
                LOGGER.debug(ids)
                LOGGER.debug(kegg_ids)
#                if len(kegg_ids) > 1:
#                    for i in range(len(kegg_ids)):
#                        mappings[databases][ids[i]] = kegg_ids[i]
#                        mappings[reverse][kegg_ids[i]] = ids[i]
#                else:
#                        mappings[databases][ids[0]] = kegg_ids[0]
#                        mappings[reverse][kegg_ids[0]] = ids[0]
        return mappings

