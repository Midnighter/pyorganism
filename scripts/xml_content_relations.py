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
import multiprocessing
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
    pool = multiprocessing.Pool()
    ask = True
    for ver in versions:
        print(ver)
        filename = os.path.join(ver, "xml_contents.pkl")
        if os.path.exists(filename):
            if ask:
                choice = raw_input("Re-parse existing XML tags? [yes|no|all|zero]")
                choice = choice.lower()
            if choice.startswith("a"):
                ask = False
            elif choice.startswith("z"):
                ask = False
                continue
            elif choice.startswith("y"):
                pass
            else:
                continue
        files = glob(os.path.join(ver, "*.xml"))
        tables = dict(pool.map(xml_reader, files))
        with open(filename, "wb") as file_h:
            pickle.dump(tables, file_h, pickle.HIGHEST_PROTOCOL)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print "Usage:\n{0} <RegulonDB base: path>".format(sys.argv[0])
        sys.exit(2)
    main(sys.argv[1])

