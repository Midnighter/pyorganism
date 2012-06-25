#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
=====================
PyOrganism Metabolism
=====================

:Authors:
    Moritz Emanuel Beber
:Date:
    2012-05-22
:Copyright:
    Copyright(c) 2012 Jacobs University of Bremen. All rights reserved.
:File:
    __init__.py
"""


from .systems import *
from .elements import *

from ..io.generic import open_file
from ..io.sbml import SBMLParser


PARSERS = {"xml": SBMLParser,
        "sbml": SBMLParser}


def read_metabolic_model(filename, frmt=False, mode="rb", encoding="utf-8", **kw_args):
    kw_args["mode"] = mode
    kw_args["encoding"] = encoding
    with  open_file(filename, **kw_args) as (file_h, ext):
        if frmt:
            ext = frmt.lower()
        if ext in PARSERS:
            parser = PARSERS[ext].get_instance()
        else:
            raise PyOrganismError("unknown metabolic system format '{0}'", ext)
        system = parser.from_string(unicode(file_h.read(-1)))
    return system

