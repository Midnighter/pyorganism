#!/usr/bin/env python
# -*- coding: utf-8 -*-


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


__all__ = ["MetabolicSystem"]


import logging
import numpy

from .. import miscellaneous as misc
from ..errors import PyOrganismError
from .elements import BasicCompartment, BasicCompound, BasicReaction


LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(misc.NullHandler())

OPTIONS = misc.OptionsManager.get_instance()


class MetabolicSystem(object):
    """
    Basically a container for reactions and compounds with some useful
    transformation functions.
    """

    def __init__(self, name="", reactions=None, compartments=None,
            compounds=None, *args, **kw_args):
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
        super(MetabolicSystem, self).__init__()
        self.name = name
        self.compartments = misc.convert(compartments, set, set())
        self.reactions = misc.convert(reactions, set, set())
        self.compounds = misc.convert(compounds, set, set())
        for rxn in self.reactions:
            self._update_compounds_compartments(rxn)
        self._transpose = None
        self._modified = True

    def __str__(self):
        return str(self.name)

    def __unicode(self):
        return unicode(self.name)

    def __eq__(self, other):
        raise NotImplementedError

    def __ne__(self, other):
        raise NotImplementedError

    def __contains__(self, element):
        if isinstance(element, BasicReaction):
            return element in self.reactions
        elif isinstance(element, BasicCompound):
            return element in self.reactions
        elif isinstance(element, BasicCompartment):
            return element in self.compartments
        else:
            raise PyOrganismError("unrecognised metabolic component '{0}'", element)

    def _update_compounds_compartments(self, reaction):
        """
        Add compounds and compartments contained within a reaction to the
        instance container.
        """
        # assumes SBMLReaction
        for cmpd in reaction.compounds():
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
        if isinstance(element, BasicReaction):
            self.reactions.add(element)
            self._update_compounds_compartments(element)
            self._modified = True
        elif isinstance(element, BasicCompound):
            self.compounds.add(element)
            self._modified = True
        elif isinstance(element, BasicCompartment):
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
        if issubclass(typeof, BasicReaction):
            self.reactions.update(elements)
            for rxn in elements:
                self._update_compounds_compartments(rxn)
            self._modified = True
        elif issubclass(typeof, BasicCompound):
            self.compounds.update(elements)
            self._modified = True
        elif issubclass(typeof, BasicCompartment):
            self.compartments.update(elements)
        else:
            raise PyOrganismError("unrecognised metabolic component type '{0}'",
                    typeof)

    def _setup_transpose(self):
        """
        Sets up a linear programming model where the transpose of the
        stoichiometric matrix is right multiplied a vector of compound masses
        and the system is expected to conform with mass conservation laws.
        """
        from ..fba import FBAModel
        if self._transpose and not self._modified:
            return
        self._transpose = FBAModel("transpose")
        # add missing reversible attribute
        for cmpd in self.compounds:
            cmpd.reversible = False
        # first add all compound masses as variables to the model
        self._transpose.add_reaction(self.compounds, lb=1.0, ub=numpy.inf)
        # constrain mass by stoichiometric coefficients
        for rxn in self.reactions:
            self._transpose.add_compound(rxn, list(rxn.compounds(True)))
#            LOGGER.debug(list(self._transpose.iter_reactions(rxn, True)))
        # objective is over all compound masses
        self._transpose.set_objective_reaction(self.compounds, 1.0)

    def verify_consistency(self, masses=False):
        """
        Verify the stoichiometric consistency of the system.

        Parameters
        ----------
        masses: bool (optional)
            If the system is consistent the minimal masses of the compounds
            should be returned.

        Returns
        -------
        bool:
            Consistent metabolic system returns True, otherwise False.
        dict:
            Optional dictionary mapping compounds to minimal masses.

        References
        ----------
        1. A. Gevorgyan, M. G Poolman, and D. A Fell, "Detection of stoichiometric
           inconsistencies in biomolecular models,"
           Bioinformatics 24, no. 19 (2008): 2245.
        """
        self._setup_transpose()
        self._transpose.fba(maximize=False)
        try:
            weights = dict(self._transpose.iter_flux())
            result = all(mass > 0.0 for mass in weights.itervalues())
        except PyOrganismError:
            LOGGER.debug(u"pssst:", exc_info=True)
            weights = dict()
            result = False
        if masses:
            return (result, weights)
        else:
            return result

    def detect_unconserved_metabolites(self):
        """
        Find those metabolites that violate the consistency of the metabolic
        system.

        Returns
        -------
        list:
            Conserved compounds in the inconsistent system.
        list:
            Unconserved compounds in the inconsistent system.

        Notes
        -----
        Before using this method make sure to verify the consistency of the
        system as this method will give wrong results for a consistent system.

        References
        ----------
        1. A. Gevorgyan, M. G Poolman, and D. A Fell, "Detection of stoichiometric
           inconsistencies in biomolecular models,"
           Bioinformatics 24, no. 19 (2008): 2245.
        """
        self._setup_transpose()
        # objective is to maximize all compound masses while they're binary
        self._transpose._make_binary(self.compounds)
        self._transpose.optimize(maximize=True)
        # sort compounds by positive or zero mass
        consistent = list()
        inconsistent = list()
        for (cmpd, value) in self._transpose.iter_flux():
            if value > 0.0:
                consistent.append(cmpd)
            else:
                inconsistent.append(cmpd)
            LOGGER.debug("%s: %f", cmpd, value)
        return (consistent, inconsistent)

    def generate_fba_model(self, name="", fluxes=False):
        """
        Generate a model fit for flux balance analysis from the metabolic
        system.
        """
        from .fba import FBAModel

        known_fluxes = list()
        objectives = list()
        factors = list()
        for rxn in self.reactions:
            if rxn.lower_bound >= 0:
                rxn.reversible = False
            if not rxn.flux_value is None:
                known_fluxes.append((rxn, rxn.flux_value))
            if rxn.objective_coefficient:
                objectives.append(rxn)
                factors.append(rxn.objective_coefficient)

        model = FBAModel(name)

        model.add_reaction(self.reactions, [list(rxn.compounds(True))\
                for rxn in self.reactions], (rxn.lower_bound\
                for rxn in self.reactions), (rxn.upper_bound\
                for rxn in self.reactions))

        model.set_objective_reaction(objectives, factors)
        if fluxes:
            return (model, dict(known_fluxes))
        else:
            return model

    def generate_network(self, disjoint_reversible=False,
            stoichiometric_coefficients=False):
        """
        Generate a network from the metabolic system.
        """
        from .networks import MetabolicNetwork
        net = MetabolicNetwork(name=self.name)
        for rxn in self.reactions:
            for cmpd in rxn.substrates:
                if stoichiometric_coefficients:
                    net.add_edge(BasicCompound(cmpd.name),
                            BasicReaction(rxn.name, rxn.reversible),
                            stoichiometry=abs(rxn.stoichiometric_coefficient(cmpd)))
                else:
                    net.add_edge(BasicCompound(cmpd.name),
                            BasicReaction(rxn.name, rxn.reversible))
            for cmpd in rxn.products:
                if stoichiometric_coefficients:
                    net.add_edge(BasicReaction(rxn.name, rxn.reversible),
                            BasicCompound(cmpd.name),
                            stoichiometry=abs(rxn.stoichiometric_coefficient(cmpd)))
                else:
                    net.add_edge(BasicReaction(rxn.name, rxn.reversible),
                            BasicCompound(cmpd.name))
            if disjoint_reversible and rxn.reversible:
                for cmpd in rxn.substrates:
                    if stoichiometric_coefficients:
                        net.add_edge(BasicReaction(
                                rxn.name + OPTIONS.reversible_suffix,
                                rxn.reversible), BasicCompound(cmpd.name),
                                stoichiometry=abs(rxn.stoichiometric_coefficient(cmpd)))
                    else:
                        net.add_edge(BasicReaction(
                                rxn.name + OPTIONS.reversible_suffix,
                                rxn.reversible), BasicCompound(cmpd.name))
                for cmpd in rxn.products:
                    if stoichiometric_coefficients:
                        net.add_edge(BasicCompound(cmpd.name),
                                BasicReaction(
                                        rxn.name + OPTIONS.reversible_suffix,
                                        rxn.reversible),
                                        stoichiometry=abs(rxn.stoichiometric_coefficient(cmpd)))
                    else:
                        net.add_edge(BasicCompound(cmpd.name),
                                BasicReaction(rxn.name, rxn.reversible))
        return net


