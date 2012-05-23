#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
=============
Model Parsers
=============

:Authors:
    Moritz Emanuel Beber
    Nikolaus Sonnenschein
:Date:
    2011-04-07
:Copyright:
    Copyright(c) 2011 Jacobs University of Bremen. All rights reserved.
:File:
    parsers.py
"""


import os
import codecs
import logging

from contextlib import contextmanager
from . import miscellaneous as misc


LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(misc.NullHandler())


def parser_warning(msg):
    LOGGER.warn("unable to parse information for:")
    LOGGER.warn(msg)

def _open_tar(path, **kw_args):
    import tarfile
    kw_args["mode"] = kw_args["mode"].strip("b")
    if isinstance(path, basestring):
        return tarfile.TarFile(name=path, mode=kw_args["mode"],
                encoding=kw_args["encoding"])
    else:
        return tarfile.TarFile(fileobj=path, mode=kw_args["mode"],
                encoding=kw_args["encoding"])

def _open_gz(path, **kw_args):
    import gzip
    if isinstance(path, basestring):
        return gzip.GzipFile(filename=path, mode=kw_args["mode"])
    else:
        return gzip.GzipFile(fileobj=path, mode=kw_args["mode"])

def _open_bz2(path, **kw_args):
    import bz2
    return bz2.BZ2File(path)

def _open_zip(path, **kw_args):
    import zipfile
    kw_args["mode"] = kw_args["mode"].strip("b")
    return zipfile.ZipFile(path, mode=kw_args["mode"])

def _open_file(path, **kw_args):
    if isinstance(path, basestring):
        return codecs.open(path, mode=kw_args["mode"],
                encoding=kw_args["encoding"])
    else:
        reader = codecs.getreader(kw_args["encoding"])
        return reader(path)

archives = {"gz": _open_gz,
        "gzip": _open_gz,
        "bz2": _open_bz2
#        "zip": _open_zip,
#        "tar": _open_tar
        }


@contextmanager
def open_file(filename, **kw_args):
    path = filename
    filename = os.path.basename(filename)
    extns = filename.split(".")
    del extns[0]
    extns.reverse()
    for ext in extns:
        ext = ext.lower()
        func = archives.get(ext, _open_file)
        path = func(path, **kw_args)
    yield (path, ext)
    if not path.closed:
        path.close()

