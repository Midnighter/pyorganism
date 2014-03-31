# -*- coding: utf-8 -*-


"""
================================
Update All Maps from One Version
================================

:Authors:
    Moritz Emanuel Beber
:Date:
    2014-01-23
:Copyright:
    Copyright(c) 2014 Jacobs University of Bremen. All rights reserved.
:File:
    update_maps.py
"""


import sys
import os
import logging

import pandas

import pyorganism
import pyorganism.regulation as pyreg # needed for gene object definition


LOGGER = logging.getLogger()
LOGGER.addHandler(logging.StreamHandler())
LOGGER.setLevel(logging.INFO)


def update_map(map_path, frame_path, update_from, version=""):
    base_path = os.path.dirname(map_path)
    if not version:
        version = os.path.basename(base_path)
    LOGGER.info("Loading genes")
    genes = pyorganism.read_pickle(os.path.join(base_path, "genes.pkl"))
    LOGGER.info("Reading data frame")
    map_name = os.path.splitext(os.path.basename(map_path))[0]
    hdf_key = "/%s" % (map_name,)
    mapping = pandas.read_hdf(frame_path, hdf_key)
#    frame_info(mapping)
    feature2gene = dict()
    for row in mapping[[version, update_from]].itertuples():
        if not row[1]:
            feature2gene[row[0]] = pyreg.Gene.get(row[2])
        elif row[1] != row[2]:
            feature2gene[row[0]] = pyreg.Gene.get(row[2])
        else:
            feature2gene[row[0]] = pyreg.Gene.get(row[1])
    pyorganism.write_pickle(feature2gene, map_path)


if __name__ == "__main__":
    if len(sys.argv) < 4 or len(sys.argv) > 5:
        LOGGER.critical("%s <feature2gene: path> <unified dataframe: path>"\
                " <update from: version> [<RegulonDB version: str>]", sys.argv[0])
        sys.exit(2)
    else:
        update_map(*sys.argv[1:])

