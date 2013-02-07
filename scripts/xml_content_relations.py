#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
=============================
Compile XML Tag Relationships
=============================

:Authors:
    Moritz Emanuel Beber
:Date:
    2012-09-06
:Copyright:
    Copyright(c) 2012 Jacobs University of Bremen. All rights reserved.
:File:
    xml_content_relations.py
"""


import sys
import os
import cPickle as pickle

from glob import glob
from bs4 import BeautifulSoup


def xml_reader(filename):
    with open(filename, "r") as file_h:
        soup = BeautifulSoup(file_h, "lxml", from_encoding="utf-8")
    name = os.path.basename(filename)
    attrs = [unicode(child.name)\
            for child in soup.row.findChildren(recursive=False)]
    return (name, attrs)

def main(db_path):
    versions = glob(os.path.join(db_path, "[0-9].[0-9]"))
    import multiprocessing
    pool = multiprocessing.Pool()
    for ver in versions:
        sys.stdout.write("\r\x1b[K")
        sys.stdout.write(ver)
        sys.stdout.flush()
        files = glob(os.path.join(ver, "*.xml"))
        tables = dict(pool.map(xml_reader, files))
        with open(os.path.join(ver, "xml_contents.pkl"), "wb") as file_h:
            pickle.dump(tables, file_h, pickle.HIGHEST_PROTOCOL)
    sys.stdout.write("\n")
    sys.stdout.flush()


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print "Usage:\n{0} <RegulonDB base: path>".format(sys.argv[0])
        sys.exit(2)
    main(sys.argv[1])

