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


__all__ = ["OptionsManager", "ProgressHandler"]


import sys
import logging

from .singletonmixin import Singleton


BREWER_SET1 = ["#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33",
        "#A65628", "#F781BF", "#8DD3C7", "#1B9E77"]


class NullHandler(logging.Handler):
    """
    A stub logging handler that ignores all logging messages. This is the
    default handler for all library loggers.
    """

    def emit(self, record):
        pass


class ProgressHandler(logging.Handler):
    """
    Write a recurring message to a stream on the same line.

    Note
    ----
    The underlying stream is not closed, as sys.stdout or sys.stderr may be used.

    Warning
    -------
    Only works properly if the message fits on a single line for the given
    stream.
    """

    def __init__(self, stream=sys.stderr, shell=False, **kw_args):
        super(ProgressHandler, self).__init__(**kw_args)
        self.stream = stream
        if hasattr(self.stream, "flush"):
            self.flush = self.stream.flush
        if shell:
            # return cursor to beginning and delete contents of line in shell
            self.head = "\r\x1b[K"
        else:
            # return cursor to beginning of line only
            self.head = "\r"

    def emit(self, record):
        """
        Emit a record.

        If a formatter is specified, it is used to format the record.
        The record is then written to the stream with a trailing newline.  If
        exception information is present, it is formatted using
        traceback.print_exception and appended to the stream.  If the stream
        has an 'encoding' attribute, it is used to determine how to do the
        output to the stream.
        """
        try:
            msg = self.format(record)
            if not logging._unicode: #if no unicode support...
                self.stream.write(self.head)
                self.stream.write(msg)
            else:
                try:
                    if (isinstance(msg, unicode) and
                        getattr(self.stream, 'encoding', None)):
                        ufs = msg.decode(self.stream.encoding)
                        try:
                            self.stream.write(self.head)
                            self.stream.write(ufs)
                        except UnicodeEncodeError:
                            #Printing to terminals sometimes fails. For example,
                            #with an encoding of 'cp1251', the above write will
                            #work if written to a stream opened or wrapped by
                            #the codecs module, but fail when writing to a
                            #terminal even when the codepage is set to cp1251.
                            #An extra encoding step seems to be needed.
                            self.stream.write(self.head)
                            self.stream.write((ufs).encode(self.stream.encoding))
                    else:
                        self.stream.write(self.head)
                        self.stream.write(msg)
                except UnicodeError:
                    self.stream.write(self.head)
                    self.stream.write(msg.encode("UTF-8"))
            self.flush()
        except StandardError:
            self.handleError(record)


class OptionsManager(Singleton):
    """
    A centralised instance to handle some common options throughout this
    library.
    """

    def __init__(self, **kw_args):
        super(OptionsManager, self).__init__(**kw_args)
        self.compound_prefix = "M_"
        self.reaction_prefix = "R_"
        self.reversible_suffix = "_Rev"
        self.exchange_reaction = "EX"
        self.compartment_suffixes = {"Cytosol": "_c", "Extra_organism": "_e",
                "Exchange": "_b", "Periplasm": "_p"}
        self.lp_solver = "gurobi"
        self.lower_bound = 0.0
        self.upper_bound = 1000.0
        self.numeric_threshold = 1E-09
        self.num_cpu = 1

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
    Convert an argument to a new type unless it is `None`.
    """
    return default if item is None else cls(item)

