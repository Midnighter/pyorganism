#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
================
Generic File I/O
================

:Authors:
    Moritz Emanuel Beber
:Date:
    2012-06-03
:Copyright:
    Copyright(c) 2012 Jacobs University of Bremen. All rights reserved.
:File:
    generic.py
"""


__all__ = ["open_file", "read_pickle", "write_pickle"]


import os
import codecs
import logging
import csv
import cPickle as pickle

from contextlib import contextmanager
from .. import miscellaneous as misc


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
    return bz2.BZ2File(path, mode=kw_args["mode"])

def _open_zip(path, **kw_args):
    import zipfile
    kw_args["mode"] = kw_args["mode"].strip("b")
    return zipfile.ZipFile(path, mode=kw_args["mode"])

def _open_file(path, **kw_args):
    if isinstance(path, basestring):
        return codecs.open(path, mode=kw_args["mode"],
                encoding=kw_args["encoding"])
    elif kw_args["encoding"]:
        if "w" in kw_args["mode"] or "a" in kw_args["mode"]:
            handle = codecs.getwriter(kw_args["encoding"])
        else:
            handle = codecs.getreader(kw_args["encoding"])
        return handle(path)
    else:
        return path

archives = {".gz": _open_gz,
        ".gzip": _open_gz,
        ".bz2": _open_bz2
#        "zip": _open_zip,
#        "tar": _open_tar
        }


@contextmanager
def open_file(filename, **kw_args):
    path = filename
    filename = os.path.basename(filename)
    while True:
        (filename, ext) = os.path.splitext(filename)
        ext = ext.lower()
        func = archives.get(ext)
        if func is None:
            break
        path = func(path, **kw_args)
    path = _open_file(path, **kw_args)
    yield (path, ext)
    if not path.closed:
        path.close()

def read_tabular(file_h, sep="\t", comment="#"):
    reader = csv.reader(file_h, delimiter=sep, quoting=csv.QUOTE_MINIMAL,
            skipinitialspace=True)
    return (row for row in reader if row and not row.startswith(comment))

def read_pickle(filename, mode="rb", **kw_args):
    kw_args["mode"] = mode.lower()
    kw_args["encoding"] = None
    with open_file(filename, **kw_args) as (file_h, ext):
        return pickle.load(file_h)

def write_pickle(instance, filename, protocol=pickle.HIGHEST_PROTOCOL,
        mode="wb", **kw_args):
    kw_args["mode"] = mode.lower()
    kw_args["encoding"] = None
    with open_file(filename, **kw_args) as (file_h, ext):
        pickle.dump(instance, file_h, protocol=protocol)

