# -*- coding: utf-8 -*-


from __future__ import (absolute_import, unicode_literals, division)


"""
============================
Metabolic Network Generators
============================

:Authors:
    Moritz Emanuel Beber
    Alexandra Mirela Grigore
:Date:
    2011-07-01
:Copyright:
    Copyright(c) 2011 Jacobs University of Bremen. All rights reserved.
:File:
    generators.py
"""


__all__ = ["random_p_mn", "random_scale_free_mn", "random_normal_scale_free"]

import logging

from builtins import (str, range)
import numpy as np

from . import elements as met
from .network import MetabolicNetwork
from .. import miscellaneous as misc


LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(misc.NullHandler())

OPTIONS = misc.OptionsManager.get_instance()


#TODO: fix seed to use new numpy rng state

def prune_network(network):
    """
    Removes stub reactions (in- or out-degree of 1 and 0 respectively) and
    assures that all other reactions consume and produce something.

    Parameters
    ----------
    network: MetabolicNetwork
        A MetabolicNetwork instance.
    """
    rand_int = np.random.random_integers
    num_rxns = len(network.reactions)
    num_cmpds = len(network.compounds)
    total = 0
    prune = list()
    for rxn in network.reactions:
        in_deg = network.in_degree(rxn)
        out_deg = network.out_degree(rxn)
        if in_deg == 0:
            if out_deg <= 1:
                prune.append(rxn)
            else:
                targets = network.successors(rxn)
                flips = rand_int(1, len(targets) - 1)
                while (flips > 0):
                    target = targets[rand_int(0, len(targets) - 1)]
                    factor = network[rxn][target]["coefficient"]
                    network.remove_edge(rxn, target)
                    network.add_edge(target, rxn, coefficient=factor)
                    LOGGER.debug("flipped direction of link %s -> %s",
                            str(rxn), str(target))
                    targets.remove(target)
                    flips -= 1
                    total += 1
        elif out_deg == 0:
            if in_deg <= 1:
                prune.append(rxn)
            else:
                targets = network.predecessors(rxn)
                flips = rand_int(1, len(targets) - 1)
                while (flips > 0):
                    target = targets[rand_int(0, len(targets) - 1)]
                    factor = network[target][rxn]["coefficient"]
                    network.remove_edge(target, rxn)
                    network.add_edge(rxn, target, coefficient=factor)
                    LOGGER.debug("flipped direction of link %s -> %s",
                            str(rxn), str(target))
                    targets.remove(target)
                    flips -= 1
                    total += 1
    for rxn in prune:
        network.remove_node(rxn)
        LOGGER.debug("removed reaction %s", str(rxn))
    prune = list()
    for cmpd in network.compounds:
        if network.degree(cmpd) == 0:
            prune.append(cmpd)
    for cmpd in prune:
        network.remove_node(cmpd)
        LOGGER.debug("removed compound %s", str(cmpd))
    LOGGER.info("%d reaction(s) and %d compound(s) removed",
            (num_rxns - len(network.reactions)),
            (num_cmpds - len(network.compounds)))
    LOGGER.info("direction of %d link(s) reversed", total)


def random_p_mn(num_compounds, num_reactions, num_reversible, p, seed=None):
    """
    Creates a bipartite graph that models metabolism according to the principles
    of an Erdos-Renyi-like random graph.

    Parameters
    ----------
    num_compounds: int
        The number of compounds (approximately) that should be in the network.
    num_reactions: int
        The number of reactions (approximately) that should be in the network.
    num_reversible: int
        The number of reactions that are reversible.
    p: float
        The probability that a link between a compound and a reaction exists.
    seed: int (optional)
        A specific seed for the random number generator for reproducible runs.

    Returns
    -------
    pymetabolism.network.networks.MetabolicNetwork

    Notes
    -----
    The numbers of compounds and reactions may be slightly less than desired
    because isolated nodes are removed from the network.
    """
    # setup
    rand_float = np.random.random_sample
    rand_int = np.random.random_integers
    cmpd_prefix = OPTIONS.compound_prefix
    rxn_prefix = OPTIONS.reaction_prefix
    if seed:
        np.random.seed(int(seed))
    num_compounds = int(num_compounds)
    num_reactions = int(num_reactions)
    num_reversible = int(num_reversible)
    p = float(p)
    network = MetabolicNetwork()
    # add compounds
    for i in range(num_compounds):
        network.add_node(met.BasicCompound("%s%d" % (cmpd_prefix, i)))
    # choose a number of reactions as reversible
    reversibles = set()
    while len(reversibles) < num_reversible:
        reversibles.add(rand_int(0, num_reactions - 1))
    for i in range(num_reactions):
        if i in reversibles:
            network.add_node(met.BasicReaction(
                    "%s%d" % (rxn_prefix, i), reversible=True))
        else:
            network.add_node(met.BasicReaction(
                "%s%d" % (rxn_prefix, i)))
    for src in network.compounds:
        for tar in network.reactions:
            if rand_float() < p:
                network.add_edge(src, tar, coefficient=0)
                LOGGER.debug("added link %s -> %s.", str(src), str(tar))
            # a conditional case here (elif not if) because we cannot determine
            # substrates and products from bidirectional edges
            elif rand_float() < p:
                network.add_edge(tar, src, coefficient=0)
                LOGGER.debug("added link %s -> %s.", str(tar), str(src))
    prune_network(network)
    return network

def random_scale_free_mn(num_compounds, num_reactions, num_reversible,
        num_rxn_tar, num_cmpd_tar, seed=None):
    """
    Uses a Barabasi-Alberts-like preferential attachment algorithm. Adopted from
    the networkx implementation.

    Parameters
    ----------
    num_compounds: int
        The number of compounds that should be in the network.
    num_reactions: int
        The number of reactions that should be in the network.
    num_reversible: int
        The number of reactions that are reversible.
    num_rxn_tar: int
        How many compounds a new reaction node links to.
    num_cmpd_tar: int
        How many reactions a new compound node links to.
    seed: int (optional)
        A specific seed for the random number generator for reproducible runs.
    """
    # setup
    rand_int = np.random.random_integers
    rand_float = np.random.random_sample
    cmpd_prefix = OPTIONS.compound_prefix
    rxn_prefix = OPTIONS.reaction_prefix
    if seed:
        np.random.seed(int(seed))
    num_compounds = int(num_compounds)
    num_reactions = int(num_reactions)
    num_reversible = int(num_reversible)
    num_rxn_tar = int(num_rxn_tar)
    num_cmpd_tar = int(num_cmpd_tar)
    network = MetabolicNetwork()
    # target nodes for reactions
    rxn_targets = []
    for i in range(num_rxn_tar):
        comp = met.BasicCompound("%s%d" % (cmpd_prefix, i))
        network.add_node(comp)
        rxn_targets.append(comp)
    # target nodes for compounds
    cmpd_targets = []
    # biased lists for preferential attachment
    repeated_cmpds = []
    repeated_rxns = []
    # choose a number of reactions as reversible
    reversibles = set()
    while len(reversibles) < num_reversible:
        reversibles.add(rand_int(0, num_reactions - 1))
    for i in range(num_cmpd_tar):
        if i in reversibles:
            rxn = met.BasicReaction("%s%d" % (rxn_prefix, i),
                    reversible=True)
        else:
            rxn = met.BasicReaction("%s%d" % (rxn_prefix, i))
        network.add_node(rxn)
        cmpd_targets.append(rxn)
        for cmpd in rxn_targets:
            if rand_float() < 0.5:
                network.add_edge(rxn, cmpd, coefficient=0)
                LOGGER.debug("added link %s -> %s", str(rxn), str(cmpd))
            else:
                network.add_edge(cmpd, rxn, coefficient=0)
                LOGGER.debug("added link %s -> %s", str(cmpd), str(rxn))
        repeated_cmpds.extend(rxn_targets)
        repeated_rxns.extend([rxn] * num_rxn_tar)
#    LOGGER.debug("Targets for compounds: %s", comp_targets)
    # current vertices being added
    current_rxn = num_cmpd_tar
    current_cmpd = num_rxn_tar
    while (current_cmpd < num_compounds or current_rxn < num_reactions):
        if current_cmpd < num_compounds:
            source = met.BasicCompound("%s%d" % (cmpd_prefix,
                    current_cmpd))
            network.add_node(source)
            for rxn in cmpd_targets:
                if rand_float() < 0.5:
                    network.add_edge(source, rxn, coefficient=0)
                    LOGGER.debug("added link %s -> %s", str(source), str(rxn))
                else:
                    network.add_edge(rxn, source, coefficient=0)
                    LOGGER.debug("added link %s -> %s", str(rxn), str(source))
            repeated_rxns.extend(cmpd_targets)
            repeated_cmpds.extend([source] * num_cmpd_tar)
            cmpd_targets = set()
            while len(cmpd_targets) < num_cmpd_tar:
                rxn = repeated_rxns[rand_int(0, len(repeated_rxns) - 1)]
                cmpd_targets.add(rxn)
            current_cmpd += 1
        if current_rxn < num_reactions:
            if current_rxn in reversibles:
                source = met.BasicReaction("%s%d" % (rxn_prefix,
                    current_rxn), reversible=True)
            else:
                source = met.BasicReaction("%s%d" % (rxn_prefix,
                    current_rxn))
            network.add_node(source)
            for cmpd in rxn_targets:
                if rand_float() < 0.5:
                    network.add_edge(source, cmpd, coefficient=0)
                    LOGGER.debug("added link %s -> %s", str(source), str(cmpd))
                else:
                    network.add_edge(cmpd, source, coefficient=0)
                    LOGGER.debug("added link %s -> %s", str(cmpd), str(source))
            repeated_cmpds.extend(rxn_targets)
            repeated_rxns.extend([source] * num_rxn_tar)
            rxn_targets = set()
            while len(rxn_targets) < num_rxn_tar:
                cmpd = repeated_cmpds[rand_int(0, len(repeated_cmpds) - 1)]
                rxn_targets.add(cmpd)
            current_rxn += 1
    prune_network(network)
    return network

def random_normal_scale_free(num_compounds, num_reactions, num_reversible,
        compound_exponent, reaction_mean, reaction_sd, seed=None):
    """
    Creates a bipartite directed graph with a normal degree distribution for
    the reaction nodes and a scale free degree distribution for the
    metabolite nodes. Adopted from the networkx implementation of
    bipartite_configuration_model.

    Parameters
    ----------
    num_compounds: int
        The number of compounds that should be in the network.
    num_reactions: int
        The number of reactions that should be in the network.
    num_reversible: int
        The number of reactions that are reversible.
    compound_exponent: float
        The exponent of the compounds' power law degree distribution.
    reaction_mean: int
        The mean of the reactions' binomial degree distribution.
    reaction_sd: float
        The standard deviation of the reactions' binomial degree distribution.
    seed: int (optional)
        A specific seed for the random number generator for reproducible runs.

    Warning
    -------
    For the given mean |mu| and standard deviation |sigma| the parameters for
    the binomial distribution are calculated as follows:

        n = |mu|^2 / (|mu| - |sigma|)
        p = (|mu| - |sigma|) / |mu|

    thus neither the mean, nor the difference between the mean and the standard
    deviation should be zero.
    """
    # setup
    rand_float = np.random.random_sample
    cmpd_prefix = OPTIONS.compound_prefix
    rxn_prefix = OPTIONS.reaction_prefix
    if seed:
        np.random.seed(int(seed))
    num_compounds = int(num_compounds)
    num_reactions = int(num_reactions)
    num_reversible = int(num_reversible)
    diff = float(reaction_mean - reaction_sd)
    if reaction_mean == 0.0 or diff == 0.0:
        raise PyOrganismError("please choose arguments that avoid"\
                " ZeroDivisionError")
    num_trials = int(round(np.power(reaction_mean, 2) / diff))
    prob = diff / reaction_mean
    compound_exponent = float(compound_exponent)
    #norm_std = int(norm_std)
    network = MetabolicNetwork()
    distr_mets = np.random.zipf(pl_exponent_compounds, num_compounds)
    distr_reacts = np.random.zipf(pl_exponent_reactions, num_reactions)
    #distr_reacts = np.random.normal(sum(distr_mets) / num_reactions, norm_std, num_reactions)
    # make the sums of the 2 degree distributions equal
    distr_mets = np.asarray(distr_mets, dtype=int)
    distr_reacts = np.asarray(distr_reacts, dtype=int)
    if sum(distr_mets) > sum(distr_reacts):
        deg_diff = sum(distr_mets) - sum(distr_reacts)
        while not deg_diff == 0:
            to_subtract_from = [n for n,i in enumerate(distr_mets) if i>2]
            h = to_subtract_from[random.randint(0, len(to_subtract_from)-1)]
            distr_mets[h] -= 1
            deg_diff -= 1
    elif sum(distr_mets) < sum(distr_reacts):
        deg_diff = sum(distr_reacts) - sum(distr_mets)
        while not deg_diff == 0:
            h = random.randint(0, num_compounds-1)
            distr_mets[h] += 1
            deg_diff -= 1

    # add metabolite nodes + build lists of degree-repeated vertices
    stubs = []
    for i in range(num_compounds):
        new_met = met.BasicCompound("%s%d" % (cmpd_prefix, i))
        network.add_node(new_met)
        stubs.extend([distr_mets[i]*[new_met]])
    astubs = [x for subseq in stubs for x in subseq]
    # add reaction  nodes + build lists of degree-repeated vertices
    stubs = []
    for i in range(num_reversible):
        new_react = met.BasicReaction("%s%d" % (rxn_prefix, i),reversible=True)
        network.add_node(new_react)
        stubs.extend([distr_reacts[i]*[new_react]])
    for i in range(num_reversible, num_reactions):
        new_react = met.BasicReaction("%s%d" % (rxn_prefix, i))
        network.add_node(new_react)
        stubs.extend([distr_reacts[i]*[new_react]])
    bstubs=[]
    bstubs=[x for subseq in stubs for x in subseq]
    # shuffle lists
    np.random.shuffle(astubs)
    np.random.shuffle(bstubs)
    # add edges
    for i in range(sum(distr_mets)):
        if rand_float() < 0.5:
            network.add_edge(astubs[i], bstubs[i], coefficient=0)
            LOGGER.debug("added link %s -> %s", str(astubs[i]), str(bstubs[i]))
        else:
            network.add_edge(bstubs[i], astubs[i], coefficient=0)
            LOGGER.debug("added link %s -> %s", str(bstubs[i]), str(astubs[i]))
            network.add_edge(u, v, coefficient=0)
            LOGGER.debug("added link %s -> %s", str(u), str(v))
    # clean up
    prune_network(network)
    return network

