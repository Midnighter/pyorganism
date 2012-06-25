#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
=============================
Reading and Writing KEGG Data
=============================

:Authors:
    Moritz Emanuel Beber
:Date:
    2012-05-11
:Copyright:
    Copyright(c) 2012 Jacobs University of Bremen. All rights reserved.
:File:
    kegg.py
"""


__all__ = ["find_organism"]


import logging
import re

from .. import miscellaneous as misc

SOAPpy = misc.load_module("SOAPpy", url="http://pywebsvcs.sourceforge.net/")


LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(misc.NullHandler())


def find_organism(self, organism, wsdl="http://soap.genome.jp/KEGG.wsdl",
        browse=10):
    """
    An interactive function that queries the KEGG Organism database and returns
    the KEGG organism identifier of the chosen result.

    Parameters
    ----------
    organism: str
        Query string for the KEGG Organism database.
    wsdl: str (optional)
        URL of the KEGG WSDL server.
    browse: int (optional)
        Maximum number of results to browse at once.

    Returns
    -------
    str:
        The three to four letter code of the chosen organism.

    Notes
    -----
    Requires SOAPpy and an active internet connection. Requires stdin.
    """
    # establish connection to DBGET server
    serv = SOAPpy.WSDL.Proxy(wsdl)
    # find appropriate organism
    # require user choice here, too many options
    choices = serv.bfind("genome " + organism)
    choices = [choice for choice in choices.split("\n") if choice]
    length = len(choices)
    start = 0
    end = min(length, browse)
    searching = True
    while start < length and searching:
        msg = [""]
        msg.append("Showing organisms %d-%d of %d, please choose an index:"\
                % (start, end - 1, length))
        msg.append("")
        for i in range(start, end):
            msg.append("[%d] %s" % (i, choices[i]))
        if end < length:
            msg.append("")
            msg.append("Type any non-integer to show the next %d organisms."\
                    % min(length - end, browse))
        msg.append("")
        try:
            selection = int(raw_input("\n".join(msg)))
        except ValueError:
            start = end
            end = min(length, end + browse)
        else:
            while True:
                try:
                    choice = choices[selection]
                except IndexError:
                    try:
                        selection = int(raw_input("Chosen index is outside"\
                                " the allowed range, try again:\n"))
                    except ValueError:
                        pass
                else:
                    searching = False
                    break
    LOGGER.info("Please be patient, this will take a few minutes.")
    pattern = re.compile(r"genome:T\d+ (\w+),")
    mobj = pattern.match(choice)
    organism = mobj.group(1)
    return organism

    def read_kegg_information(self, organism,
            wsdl="http://soap.genome.jp/KEGG.wsdl", num_threads=20):
        """
        A threaded method that extracts reactions information from KEGG
        pathways.

        Compound-reaction pairs are added as links to the network. A
        link attribute "rpair" contains the reaction pair information "main",
        "trans", or "leave".

        Parameters
        ----------
        organism: str
            KEGG Organism identifier consisting of 3-4 lower case letters.
        wsdl: str (optional)
            URL of the KEGG WSDL server.
        num_threads: int (optional)
            The number of desired simultaneous connections to the KEGG WSDL
            server.

        Notes
        -----
        Requires SOAPpy and an active internet connection.
        """
        from .wsdl import ThreadedWSDLFetcher
        from Queue import Queue
        # establish connection to DBGET server
        serv = SOAPpy.WSDL.Proxy(wsdl)
        pathways = serv.list_pathways(organism)
        LOGGER.info("KEGG contains {0:d} pathways for the organism '{1}'.",
                len(pathways), organism)
        # use a threaded approach to server querying
        tasks = Queue()
        for i in range(num_threads):
            thrd = ThreadedWSDLFetcher(tasks, wsdl)
            thrd.start()
        reactions = list()
        for path in pathways:
            tasks.put(("get_reactions_by_pathway", path.entry_id, reactions))
        tasks.join()
        reactions = set([rxn for objs in reactions for rxn in objs])
        LOGGER.info("The pathways contain %d unique reactions", len(reactions))
        descriptions = list()
        for rxn in reactions:
            tasks.put(("bget", rxn, descriptions))
        tasks.join()
        return descriptions

