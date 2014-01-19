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
from lxml import etree


LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.INFO)
LOGGER.addHandler(logging.StreamHandler())


def xml_reader(filename):
    """
    A method using iterparse as above would be preferable, since we just want to
    collect the first few tags. Unfortunately, so far iterparse does not work
    with html (aka broken xml).
    """
    name = os.path.basename(filename)
    with open(filename, "rb") as file_h:
        if etree.LXML_VERSION < (3, 3):
            parser = etree.HTMLParser(encoding="latin1")
            tree = etree.parse(file_h, parser)
            row_it = tree.iter(tag="row")
            element = next(row_it)
            attrs = [unicode(child.tag) for child in element.iterchildren()]
        else:
            row_it = etree.iterparse(file_h, tag="row", html=True)
            (event, element) = next(row_it)
            attrs = [unicode(child.tag) for child in element.iterchildren()]
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

