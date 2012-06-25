#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
==================================
RegulonDB Parsing and Web Services
==================================

:Authors:
    Moritz Emanuel Beber
:Date:
    2012-06-03
:Copyright:
    Copyright(c) 2012 Jacobs University of Bremen. All rights reserved.
:File:
    regulondb.py
"""


import logging
import re

from .. import miscellaneous as misc
from ..regulation.elements import Gene, TranscriptionalRegulator


LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(misc.NullHandler())


def parse_regulondb_entry(eck12, dscrptn):
    dscrptn = dscrptn.split("\n")
    gene = False
    product = False
    info = dict()
    prod_info = dict()
    blattner = re.compile(r"b\d+")
    location = re.compile(r"\d+")
    phantom = 0
    for line in dscrptn:
        if line.startswith("GENE"):
            gene = True
        elif line.startswith("PRODUCT"):
            gene = False
            product = True
        elif gene and line.startswith("\t"):
            line = line.strip()
            tmp = line.split(":")
            if tmp[0].startswith("Name"):
                name = tmp[1].strip()
                if name == "Phantom Gene":
                    phantom += 1
                    info["name"] = "%s %d" % (name, phantom)
                else:
                    info["name"] = name
            elif tmp[0].startswith("Synonym"):
                info["synonyms"] = list()
                for name in tmp[1].split(","):
                    name = name.strip()
                    if blattner.match(name):
                        info["bnumber"] = name
                    else:
                        info["synonyms"].append(name)
            elif tmp[0].startswith("Map"):
                info["position"] = location.findall(tmp[1])
                info["position"] = map(int, info["position"])
                info["position"] = tuple(info["position"])
            elif tmp[0].startswith("Strand"):
                info["strand"] = tmp[1].strip()
        elif product and line.startswith("\t"):
            line = line.strip()
            tmp = line.split(":")
            if tmp[0].startswith("Name"):
                prod_info["name"] = tmp[1].strip()
            elif tmp[0].startswith("Synonym"):
                prod_info["synonyms"] = list()
                for name in tmp[1].split(","):
                    name = name.strip()
                    prod_info["synonyms"].append(name)
            elif tmp[0].startswith("GO"):
                prod_info["go"] = int(tmp[2])
        else:
            gene = False
            product = False
            break
    tf = TranscriptionalRegulator(**prod_info)
    return Gene(eck12=eck12, product=tf, **info)

def read_gene_information(eck12_nums, num_threads=20,
        wsdl="http://regulondb.ccg.unam.mx/webservices/Gene.jws?wsdl"):
    """
    A threaded method that extracts gene information from RegulonDB.

    Parameters
    ----------
    eck12_nums: list
        A collection of ECK12 identifiers used by RegulonDB whose information
        should be retrieved.
    wsdl: str (optional)
        URL of the RegulonDB WSDL server.
    num_threads: int
        The number of desired simultaneous connections to the WSDL server.

    Notes
    -----
    Requires SOAPpy and an active internet connection.
    """

    from .wsdl import ThreadedWSDLFetcher
    from Queue import Queue
    # use a threaded approach to server querying
    tasks = Queue()
    for i in range(num_threads):
        thrd = ThreadedWSDLFetcher(tasks, wsdl)
        thrd.start()
    # get all gene descriptions from RegulonDB
    descriptions = list()
    for gene in eck12_nums:
        tasks.put(("getGene", gene, descriptions))
    tasks.join()
    return descriptions

