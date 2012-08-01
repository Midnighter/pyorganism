#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
=====================
PyOrganism Metabolism
=====================

:Authors:
    Moritz Emanuel Beber
:Date:
    2012-05-22
:Copyright:
    Copyright(c) 2012 Jacobs University of Bremen. All rights reserved.
:File:
    __init__.py
"""


import logging

from .. import miscellaneous as misc
from ..errors import PyOrganismError
from .systems import MetabolicSystem
from .networks import MetabolicNetwork
from ..io.generic import open_file
from ..io.sbml import SBMLParser
from .elements import *
from ..io.bigg import BiGGCompoundParser
from ..io.metrxn import *


LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(misc.NullHandler())

OPTIONS = misc.OptionsManager.get_instance()

PARSERS = {".xml": SBMLParser,
        ".sbml": SBMLParser}


def read_metabolic_model(filename, frmt=None, mode="rb", encoding="utf-8", **kw_args):
    kw_args["mode"] = mode
    kw_args["encoding"] = encoding
    with  open_file(filename, **kw_args) as (file_h, ext):
        if not frmt is None:
            ext = frmt.lower()
        if ext in PARSERS:
            parser = PARSERS[ext].get_instance()
        else:
            raise PyOrganismError("unknown metabolic system format '{0}'", ext)
        system = parser.from_string(str(file_h.read(-1)))
    return system

def read_edgelist(filename, delimiter=None, comments="#", mode="rb",
        encoding="utf-8", **kw_args):
    """
    """
    def build_node(name):
        if name.startswith(OPTIONS.compound_prefix):
            compound = BasicCompound(name[len(OPTIONS.compound_prefix):])
            return compound
        elif name.startswith(OPTIONS.reaction_prefix):
            if name.endswith(OPTIONS.reversible_suffix):
                reaction = BasicReaction(name[len(OPTIONS.reaction_prefix):
                        -len(OPTIONS.reversible_suffix)])
                reaction.reversible = True
            else:
                reaction = BasicReaction(name[len(OPTIONS.reaction_prefix):])
            return reaction
        else:
            raise PyOrganismError("unidentified metabolic object '{0}'".format(name))

    # read the file contents
    kw_args["mode"] = mode
    kw_args["encoding"] = encoding
    with  open_file(filename, **kw_args) as (file_h, ext):
        lines = [line.strip() for line in file_h]
    # parse the edgelist into a simple metabolic network
    net = MetabolicNetwork(name=filename)
    for line in lines:
        line = line.strip()
        if line.startswith(comments) or line == "":
            continue
        tmp = line.split(delimiter)
        u = build_node(tmp[0])
        if isinstance(u, BasicReaction) and\
                tmp[0].endswith(OPTIONS.reversible_suffix):
            continue
        v = build_node(tmp[1])
        if isinstance(v, BasicReaction) and\
                tmp[1].endswith(OPTIONS.reversible_suffix):
            continue
        net.add_edge(u, v)
    return net

def write_edgelist(network, filename, distinct=True, delimiter="\t",
        comments="#", mode="wb", encoding="utf-8", **kw_args):
    """
    """
    # assemble lines
    for rxn in network.reactions:
        rxn_name = OPTIONS.reaction_prefix + rxn.name
        if rxn.reversible:
            if distinct:
                rev_name = "{0}{1}{2}".format(OPTIONS.reaction_prefix, rxn.name,
                        OPTIONS.reversible_suffix)
            else:
                rev_name = rxn_name
            for cmpd in self.successors_iter(rxn):
                lines.append("{0}{1}{2}\n".format(rxn_name, delimiter,
                        OPTIONS.compound_prefix + cmpd.name))
                lines.append("{0}{1}{2}\n".format(OPTIONS.compound_prefix + cmpd.name,
                        delimiter, rev_name))
            for cmpd in self.predecessors_iter(rxn):
                lines.append("{0}{1}{2}\n".format(OPTIONS.compound_prefix + cmpd.name,
                        delimiter, rxn_name))
                lines.append("{0}{1}{2}\n".format(rev_name, delimiter,
                        OPTIONS.compound_prefix + cmpd.name))
        else:
            for cmpd in self.successors_iter(rxn):
                lines.append("{0}{1}{2}\n".format(rxn_name, delimiter,
                        OPTIONS.compound_prefix + cmpd.name))
            for cmpd in self.predecessors_iter(rxn):
                lines.append("{0}{1}{2}\n".format(OPTIONS.compound_prefix + cmpd.name,
                        delimiter, rxn_name))
    # write to file
    kw_args["mode"] = mode
    kw_args["encoding"] = encoding
    with  open_file(filename, **kw_args) as (file_h, ext):
        file_handle.writelines(lines)

    def read_kegg_reactions(descriptions):
        pattern = re.compile(r"\S")
        for info in descriptions:
            info = info.split("\n")
            name = info[0].split()[1]
            begin = -1
            stop = 0
            for (i, line) in enumerate(info):
                if line.startswith("RPAIR"):
                    begin = i
                    continue
                if begin > -1 and pattern.match(line):
                    stop = i
                    break
            if begin < 0:
                LOGGER.warn("No reaction pair information for '%s'.", name)
                continue
            pairs = [info[begin].split()[1:]]
            for i in range(begin + 1, stop):
                pairs.append(info[i].split())
            LOGGER.debug(str(pairs))
            reac = BasicReaction(name, reversible=True)
            for line in pairs:
                (u, v) = line[1].split("_")
                self.add_edge(BasicCompound(u), reac, rpair=line[2])
                self.add_edge(reac, BasicCompound(v), rpair=line[2])

