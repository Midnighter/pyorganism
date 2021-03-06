# -*- coding: utf-8 -*-


from __future__ import (absolute_import, unicode_literals)


"""
=================
Metabolic Systems
=================

:Authors:
    Moritz Emanuel Beber
    Nikolaus Sonnenschein
:Date:
    2011-04-07
:Copyright:
    Copyright(c) 2011 Jacobs University of Bremen. All rights reserved.
:File:
    systems.py
"""


__all__ = ["MetabolicSystem", "read_metabolic_model", "generate_fba_model"]


import logging

from builtins import (str, dict)
from future.utils import python_2_unicode_compatible

from .. import miscellaneous as misc
from ..errors import PyOrganismError
from ..io.generic import open_file
from ..io.sbml import SBMLParser
from . import elements as pymet


LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(misc.NullHandler())

OPTIONS = misc.OptionsManager.get_instance()

PARSERS = {".xml": SBMLParser,
        ".sbml": SBMLParser}


@python_2_unicode_compatible
class MetabolicSystem(object):
    """
    Basically a container for reactions and compounds with some useful
    transformation functions.
    """

    def __init__(self, name="", reactions=None, compartments=None,
            compounds=None, **kw_args):
        """
        Parameters
        ----------
        name: str (optional)
            A string uniquely identifying the `MetabolicSystem` instance.
        reactions: iterable (optional)
            An iterable of `BasicReaction` instances. `BasicCompound`s and
            `BasicCompartments` within those reactions are automatically added
            to the system.
        compartments: iterable (optional)
            An iterable of `BasicCompartment` instances.
        compounds: iterable (optional)
            Additional compounds not contained in the reactions that should be
            added to the system.
        """
        super(MetabolicSystem, self).__init__(**kw_args)
        self.name = str(name)
        self.compartments = misc.convert(compartments, set, set())
        self.reactions = misc.convert(reactions, set, set())
        self.compounds = misc.convert(compounds, set, set())
        for rxn in self.reactions:
            self._update_compounds_compartments(rxn)
        self._transpose = None
        self._modified = True

    def __str__(self):
        return self.name

    def __eq__(self, other):
        return NotImplemented

    def __ne__(self, other):
        return NotImplemented

    def __contains__(self, element):
        if isinstance(element, pymet.BasicReaction):
            return element in self.reactions
        elif isinstance(element, pymet.BasicCompound):
            return element in self.reactions
        elif isinstance(element, pymet.BasicCompartment):
            return element in self.compartments
        else:
            raise PyOrganismError("unrecognised metabolic component '{0}'", element)

    def _update_compounds_compartments(self, reaction):
        """
        Add compounds and compartments contained within a reaction to the
        instance container.
        """
        # assumes SBMLReaction
        for cmpd in reaction.compounds_iter():
            self.compounds.add(cmpd)
            if hasattr(cmpd, "compartment"):
                self.compartments.add(cmpd.compartment)

    def add(self, element):
        """
        Adds a single compartment, reaction, or compound to the metabolic
        system.

        Parameters
        ----------
        element:
            A child instance of BasicMetabolicComponent.
        """
        if isinstance(element, pymet.BasicReaction):
            self.reactions.add(element)
            self._update_compounds_compartments(element)
            self._modified = True
        elif isinstance(element, pymet.BasicCompound):
            self.compounds.add(element)
            self._modified = True
        elif isinstance(element, pymet.BasicCompartment):
            self.compartments.add(element)
        else:
            raise PyOrganismError("unrecognised metabolic component type '{0}'",
                    type(element))

    def update(self, elements, typeof):
        """
        Adds a compartments, reactions, or compounds to the metabolic
        system.

        Parameters
        ----------
        elements: iterable
            Iterable of BasicMetabolicComponents.
        typeof: class
            The specific type of the elements.

        Warning
        -------
        If elements contains different metabolic components you risk messing up the
        class' internal structure. In that case rather use the add function for
        each component. This is a convenience function only.
        """
        elements = set(elements)
        if issubclass(typeof, pymet.BasicReaction):
            self.reactions.update(elements)
            for rxn in elements:
                self._update_compounds_compartments(rxn)
            self._modified = True
        elif issubclass(typeof, pymet.BasicCompound):
            self.compounds.update(elements)
            self._modified = True
        elif issubclass(typeof, pymet.BasicCompartment):
            self.compartments.update(elements)
        else:
            raise PyOrganismError("unrecognised metabolic component type '{0}'",
                    typeof)

    def remove_compartment(self, compartment):
        rm = set()
        if hasattr(compartment, "__iter__"):
            for rxn in self.reactions:
                if any(cmpd in cmprtmnt for cmpd in rxn.compounds_iter() for
                        cmprtmnt in compartment):
                    rm.add(rxn)
            self.reactions.difference_update(rm)
            rm = set([cmpd for cmprtmnt in compartment for cmpd in cmprtmnt])
            self.compounds.difference_update(rm)
            self.compartments.difference_update(set(compartment))
        else:
            for rxn in self.reactions:
                if any(cmpd in compartment for cmpd in rxn.compounds_iter()):
                    rm.add(rxn)
            self.reactions.difference_update(rm)
            rm = set([cmpd for cmpd in compartment])
            self.compounds.difference_update(rm)
            self.compartments.remove(compartment)

    def decompartmentalize(self):
        for rxn in self.reactions:
            members = dict()
            for (cmpd, factor) in rxn.substrates.items():
                if isinstance(cmpd, pymet.BasicCompartmentCompound):
                    members[cmpd.compound] = factor
                else:
                    members[cmpd] = factor
            rxn.substrates = members
            members = dict()
            for (cmpd, factor) in rxn.products.items():
                if isinstance(cmpd, pymet.BasicCompartmentCompound):
                    members[cmpd.compound] = factor
                else:
                    members[cmpd] = factor
        self.compartments = set()

    def generate_network(self, disjoint_reversible=False,
            stoichiometric_coefficients=False):
        """
        Generate a network from the metabolic system.
        """
        from .networks import MetabolicNetwork
        net = MetabolicNetwork(name=self.name, compartments=self.compartments)
        for cmpd in self.compounds:
            net.add_node(cmpd)
        for rxn in self.reactions:
            net.add_node(rxn)
            for cmpd in rxn.substrates:
                if stoichiometric_coefficients:
                    net.add_edge(cmpd, rxn,
                            stoichiometry=rxn.stoichiometric_coefficient(cmpd))
#                    if rxn.reversible:
#                        net.add_edge(rxn, cmpd,
#                                stoichiometry=rxn.stoichiometric_coefficient(cmpd))
                else:
                    net.add_edge(cmpd, rxn)
#                    if rxn.reversible:
#                        net.add_edge(rxn, cmpd)
            for cmpd in rxn.products:
                if stoichiometric_coefficients:
                    net.add_edge(rxn, cmpd,
                            stoichiometry=rxn.stoichiometric_coefficient(cmpd))
#                    if rxn.reversible:
#                        net.add_edge(cmpd, rxn,
#                                stoichiometry=rxn.stoichiometric_coefficient(cmpd))
                else:
                    net.add_edge(rxn, cmpd)
#                    if rxn.reversible:
#                        net.add_edge(cmpd, rxn)
#            if disjoint_reversible and rxn.reversible:
#                for cmpd in rxn.substrates:
#                    if stoichiometric_coefficients:
#                        net.add_edge(pymet.BasicReaction(
#                                rxn.name + OPTIONS.reversible_suffix,
#                                rxn.reversible), cmpd,
#                                stoichiometry=abs(rxn.stoichiometric_coefficient(cmpd)))
#                    else:
#                        net.add_edge(pymet.BasicReaction(
#                                rxn.name + OPTIONS.reversible_suffix,
#                                rxn.reversible), cmpd)
#                for cmpd in rxn.products:
#                    if stoichiometric_coefficients:
#                        net.add_edge(cmpd,
#                                pymet.BasicReaction(
#                                        rxn.name + OPTIONS.reversible_suffix,
#                                        rxn.reversible),
#                                        stoichiometry=abs(rxn.stoichiometric_coefficient(cmpd)))
#                    else:
#                        net.add_edge(cmpd, rxn)
        return net

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

def generate_fba_model(metabolism, name="", fluxes=False):
    """
    Generate a model fit for flux balance analysis from the metabolic
    system.
    """
    from .fba import FBAModel

    known_fluxes = list()
    objectives = list()
    factors = list()
    for rxn in metabolism.reactions:
        if not rxn.flux_value is None:
            known_fluxes.append((rxn, rxn.flux_value))
        if rxn.objective_coefficient:
            objectives.append(rxn)
            factors.append(rxn.objective_coefficient)
    model = FBAModel(name)
    model.add_reaction(metabolism.reactions, [list(rxn.compounds_iter(True))\
            for rxn in metabolism.reactions], (rxn.lower_bound\
            for rxn in metabolism.reactions), (rxn.upper_bound\
            for rxn in metabolism.reactions))
    model.set_objective_reaction(objectives, factors)
    exchange = pymet.SBMLCompartment["EX"]
    for cmpd in exchange.iter_compartmentalized():
        model.add_drain(cmpd, lb=OPTIONS.lower_bound, ub=OPTIONS.upper_bound)
        model.add_source(cmpd, lb=OPTIONS.lower_bound, ub=OPTIONS.upper_bound)
    if fluxes:
        return (model, dict(known_fluxes))
    else:
        return model

