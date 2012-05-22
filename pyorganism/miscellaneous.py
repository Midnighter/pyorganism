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


__all__ = ["OptionsManager"]


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
        Singleton.__init__(self, *args, **kw_args)
        self.num_proc = 1

def load_module(module, name=False, url=False):
    if not sys.modules.has_key(module):
        try:
            __import__(module)
        except ImportError:
            if not name:
                name = module
            msg = list()
            msg.append(u"%s is required for this functionality." % name)
            msg.append(u"Please specify a different external dependency in the")
            msg.append(u"options or double-check your installation and necessary")
            msg.append(u"Python bindings.")
            if url:
                msg.append(u"Please see %s for detailed installation" % url)
                msg.append(u"instructions.")
            raise ImportError(" ".join(msg))
    return sys.modules[module]


