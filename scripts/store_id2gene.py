# -*- coding: utf-8 -*-


"""
===============================
Unified Storage of Gene ID Maps
===============================

:Authors:
    Moritz Emanuel Beber
:Date:
    2014-01-23
:Copyright:
    Copyright(c) 2014 Jacobs University of Bremen. All rights reserved.
:File:
    store_id2gene.py
"""


import sys
import os
import logging
import warnings

import pandas

import pyorganism
import pyorganism.regulation as pyreg # needed for gene object definition


LOGGER = logging.getLogger()
LOGGER.addHandler(logging.StreamHandler())
LOGGER.setLevel(logging.INFO)


def create_dataframe(feature2gene, version):
    features = feature2gene.keys()
    unique_ids = [str(feature2gene[item].unique_id) if feature2gene[item] else ""\
            for item in features]
    df = pandas.DataFrame(data=unique_ids, index=features, columns=(version,))
    df.sort(inplace=True)
    return df

def frame_info(df):
    for col in df.columns:
        LOGGER.info("Version %s:", col)
        num = (df[col] == "").sum()
        LOGGER.info("    %d unmapped features", num)

def unequal_row(row):
    return not all(x == row[0] for x in row[1:])

def missing_entry(row):
    return any(x == "" for x in row)

def inconsistencies(df):
    # using itertuples is faster than using df.apply(..., axis=1, raw=True)
    # chaining the functions but does not return a Series object
    faulty_indeces = [unequal_row(row) or missing_entry(row) for row in
            df.itertuples(index=False)]
    return faulty_indeces

def store_joint_ids(map_path, frame_path, version=""):
    if not version:
        version = os.path.basename(map_path)
    base_path = os.path.dirname(map_path)
    if not version:
        version = os.path.basename(base_path)
    LOGGER.info("{0:*^78s}".format(version))
    # load into memory
    LOGGER.info("Loading genes")
    genes = pyorganism.read_pickle(os.path.join(base_path, "genes.pkl"))
    LOGGER.info("Loading feature map")
    id2gene = pyorganism.read_pickle(map_path)
    LOGGER.info("Creating data frame")
    df = create_dataframe(id2gene, version)
    map_name = os.path.splitext(os.path.basename(map_path))[0]
    hdf_key = "/%s" % (map_name,)
    LOGGER.info("Assembling joint frame")
    if os.path.exists(frame_path):
        try:
            mapping = pandas.read_hdf(frame_path, hdf_key)
            mapping[version] = pandas.Series(df[version], index=mapping.index)
        except KeyError:
            mapping = df
    else:
        mapping = df
    frame_info(mapping)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", pandas.io.pytables.PerformanceWarning)
        mapping.to_hdf(frame_path, hdf_key, format="table")



if __name__ == "__main__":
    if len(sys.argv) < 3 or len(sys.argv) > 4:
        LOGGER.critical("%s <feature2gene: path> <unified dataframe: path>"\
                " [<RegulonDB version: str>]", sys.argv[0])
        sys.exit(2)
    else:
        store_joint_ids(*sys.argv[1:])

