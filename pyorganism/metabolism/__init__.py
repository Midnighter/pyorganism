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


import warnings

from .elements import *
from .systems import *
from .networks import *
try:
    from .fba import *
    from .nullmodels import *
except ImportError:
    warnings.warn("It seems you have no supported linear solver installed\
            (Gurobi): FBA is disabled.")
from .generators import *
try:
    from ..io.sbml import *
except ImportError:
    warnings.warn("You do not have the required libsbml package: SBML file io\
            is disabled.")
from ..io.bigg import *
from ..io.metrxn import *
try:
    from ..io.kegg import *
except ImportError:
    warnings.warn("You do not have the required SOAPpy package: KEGG interface\
            is disabled.")

def clear():
    """
    Clear all memories of elements' classes.
    """
    pass
