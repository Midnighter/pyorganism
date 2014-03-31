# -*- coding: utf-8 -*-


import sys
import os
import logging

import pandas

import pyorganism
import pyorganism.regulation as pyreg


LOGGER = logging.getLogger()
LOGGER.addHandler(logging.StreamHandler())
LOGGER.setLevel(logging.INFO)


def construct_name2gene(objects_path, map_path):
    version = os.path.basename(objects_path)
    if not version:
        version = os.path.basename(os.path.dirname(objects_path))
    # load objects
    pyorganism.read_pickle(os.path.join(objects_path, "genes.pkl"))
    mapping = pandas.read_hdf(map_path, "/corrected")
    name2gene = dict((name, pyreg.Gene.get(eck12)) for (name, eck12) in
            mapping[version].iteritems())
    num = sum(1 for gene in name2gene.itervalues() if gene)
    LOGGER.info("%d names were mapped onto %d genes (%3.2f%%).", len(name2gene),
            num, 100.0 * num / len(name2gene))
    pyorganism.write_pickle(name2gene, os.path.join(objects_path, "name2gene.pkl"))

if __name__ == "__main__":
    if len(sys.argv) != 3:
        LOGGER.critical("%s <RegulonDB objects path> <map path>", sys.argv[0])
        sys.exit(2)
    else:
        construct_name2gene(sys.argv[1], sys.argv[2])

