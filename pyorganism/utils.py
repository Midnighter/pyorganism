# -*- coding: utf-8 -*-


from __future__ import (absolute_import, unicode_literals)


"""
==================
PyOrganism Utility
==================

:Authors:
    Moritz Emanuel Beber
:Date:
    2014-11-01
:Copyright:
    Copyright(c) 2014 Jacobs University of Bremen. All rights reserved.
:File:
    utils.py
"""


__all__ = ["version_from_path"]


import logging
from os.path import (basename, dirname)

from . import miscellaneous as misc


LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(misc.NullHandler())


def version_from_path(path):
    ver = basename(path)
    if not ver:
        ver = basename(dirname(path))
    return ver

