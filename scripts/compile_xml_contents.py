#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
========================
Compile XML Tag Contents
========================

:Authors:
    Moritz Emanuel Beber
:Date:
    2012-09-06
:Copyright:
    Copyright(c) 2012 Jacobs University of Bremen. All rights reserved.
:File:
    compile_xml_contents.py
"""


import logging
import sys
import os
import multiprocessing
import cPickle as pickle

from glob import glob
from bs4 import BeautifulSoup


LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.INFO)
LOGGER.addHandler(logging.StreamHandler())


def xml_reader(filename):
    with open(filename, "r") as file_h:
        soup = BeautifulSoup(file_h, "lxml", from_encoding="utf-8")
    name = os.path.basename(filename)
    attrs = [unicode(child.name)\
            for child in soup.row.findChildren(recursive=False)]
    return (name, attrs)

def main(db_path, out_dir):
    versions = glob(os.path.join(db_path, "[0-9].[0-9]"))
    pool = multiprocessing.Pool()
    for ver in versions:
        ver_name = os.path.basename(ver)
        LOGGER.info(ver_name)
        out_sub = os.path.join(out_dir, ver_name)
        out_name = os.path.join(out_sub, "xml_contents.pkl")
        if os.path.exists(out_name):
            continue
        if not os.path.exists(out_sub):
            os.makedirs(out_sub)
        files = glob(os.path.join(ver, "*.xml"))
        tables = dict(pool.map(xml_reader, files))
        with open(out_name, "wb") as file_h:
            pickle.dump(tables, file_h, pickle.HIGHEST_PROTOCOL)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print "Usage:\n{0} <RegulonDB base: path> <output: path>".format(sys.argv[0])
        sys.exit(2)
    main(sys.argv[1], sys.argv[2])

