#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
===================
Library Miscellanea
===================

:Authors:
    Moritz Emanuel Beber
:Date:
    2012-05-22
:Copyright:
    Copyright(c) 2012 Jacobs University of Bremen. All rights reserved.
:File:
    miscellaneous.py
"""


__all__ = ["OptionsManager", "convert"]


import sys
import logging

from .singletonmixin import Singleton


class NullHandler(logging.Handler):
    """
    A stub logging handler that ignores all logging messages. This is the
    default handler for all library loggers.
    """

    def emit(self, record):
        pass


class OptionsManager(Singleton):
    """
    A centralised instance to handle some common options throughout this
    library.
    """

    def __init__(self, *args, **kw_args):
        super(OptionsManager, self).__init__(*args, **kw_args)
        self.num_proc = 1

def load_module(module, name=False, url=False):
    try:
        mod = sys.modules.get(module, __import__(module))
    except ImportError:
        if not name:
            name = module
        msg = list()
        msg.append(u"{0} is required for this functionality.".format(name))
        msg.append(u"Please specify a different external dependency in the")
        msg.append(u"options or double-check your installation and necessary")
        msg.append(u"Python bindings.")
        if url:
            msg.append(u"Please see {0} for detailed installation".format(url))
            msg.append(u"instructions.")
        raise ImportError(" ".join(msg))
    return mod


def convert(item, cls, default=None):
    """
    Convert the item to the specified type unless item is `None`.
    """
    return default if item is None else cls(item)

