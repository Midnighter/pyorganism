#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
===========================
SBML Metabolic Model Parser
===========================

:Authors:
    Moritz Emanuel Beber
    Nikolaus Sonnenschein
:Date:
    2011-04-07
:Copyright:
    Copyright(c) 2011 Jacobs University of Bremen. All rights reserved.
:File:
    sbml.py
"""


__all__ = ["SBMLParser"]


import logging

from ..singletonmixin import Singleton
from .. import miscellaneous as misc
from ..metabolism import elements as pymet
from ..metabolism.systems import MetabolicSystem


LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(misc.NullHandler())

OPTIONS = misc.OptionsManager.get_instance()


class SBMLParser(Singleton):
    """
    A class implementing methods for parsing a SBML model
    """

    _sbml = False

    def __init__(self, **kw_args):
        if not self.__class__._sbml:
            self.__class__._sbml = misc.load_module("libsbml", "SBML",
                    "http://sbml.org/Software/libSBML")
        super(SBMLParser, self).__init__(**kw_args)

    def from_string(self, xml):
        """
        Parse a document in SBML format.
        """
        document = self._sbml.readSBMLFromString(xml)
        if document.getNumErrors() > 0:
            LOGGER.warn("reading the SBML document produced some errors")
        model = document.getModel()
        # parse compartments
        compartments = [self._parse_compartment(comp) for comp in
                model.getListOfCompartments()]
        LOGGER.debug("approx. %d compartments", len(compartments))
        # parse compounds
        compounds = [self._parse_species(cmpd) for cmpd in
                model.getListOfSpecies()]
        LOGGER.debug("%d compounds", len(compounds))
        reactions = [self._parse_reaction(rxn, model) for rxn in
                model.getListOfReactions()]
        LOGGER.debug("%d reactions", len(reactions))
        return MetabolicSystem(compartments=compartments,
                reactions=reactions, compounds=compounds)

    def _parse_compartment(self, compartment):
        suffix = ""
        for (suff, name) in OPTIONS.compartments.iteritems():
            if name == compartment.getId():
                suffix = suff
        return pymet.SBMLCompartment(unique_id=compartment.getId(),
                name=compartment.getName(),
                outside=compartment.getOutside(),
                constant=compartment.getConstant(),
                spatial_dimensions=compartment.getSpatialDimensions(),
                size=compartment.getSize(), units=compartment.getUnits(),
                suffix=suffix)

    def _strip_species_id(self, name):
        identifier = name
        if identifier.startswith(OPTIONS.compound_prefix):
            identifier = identifier[len(OPTIONS.compound_prefix):]
        compartment = None
        for suffix in OPTIONS.compartments.iterkeys():
            if identifier.endswith(suffix):
                identifier = identifier[:-len(suffix)]
                compartment = pymet.SBMLCompartment(
                        unique_id=OPTIONS.compartments[suffix], suffix=suffix)
                break
        return (identifier, compartment)

    def _parse_species(self, compound):
        """
        Able to parse entries from getListOfSpecies

        @todo: Check for meta information and parse if available
        """
        (identifier, comp) = self._strip_species_id(compound.getId())
        if comp is None:
            comp = pymet.SBMLCompartment(compound.getCompartment())
        name = compound.getName()
        cmpd = pymet.SBMLCompound(unique_id=identifier, name=name,
                charge=compound.getCharge())
        if not comp:
            return cmpd
        else:
            return pymet.SBMLCompartmentCompound(unique_id=identifier + comp.suffix,
                    compound=cmpd, compartment=comp)

    def _strip_reaction_id(self, name):
        identifier = name
        if identifier.startswith(OPTIONS.reaction_prefix):
            identifier = identifier[len(OPTIONS.reaction_prefix):]
        if identifier.endswith(OPTIONS.reversible_suffix):
            identifier = identifier[:-len(OPTIONS.reversible_suffix)]
        return identifier

    def _parse_reaction(self, reaction, model, note_sep=":"):
        """Able to parse entries from getListOfReactions"""
        identifier = self._strip_reaction_id(reaction.getId())
        name = reaction.getName()
        # parse additional reaction parameters
        params = dict()
        for param in reaction.getKineticLaw().getListOfParameters():
            params[param.getId().lower()] = param.getValue()
        # substrates' stoichiometry
        substrates = dict((self._parse_species(model.getSpecies(
            elem.getSpecies())),
            abs(elem.getStoichiometry())) for elem in
            reaction.getListOfReactants())
        # products' stoichiometry
        products = dict((self._parse_species(model.getSpecies(
            elem.getSpecies())),
            abs(elem.getStoichiometry())) for elem in
            reaction.getListOfProducts())
        # other information contained in notes
        info = dict()
        notes = reaction.getNotes()
        for i in range(notes.getNumChildren()):
            node = notes.getChild(i)
            for j in range(node.getNumChildren()):
                item = node.getChild(j).toString().split(note_sep, 1)
                if len(item) == 2:
                    key = item[0].strip().lower().replace(" ", "_")
                    value = item[1].strip()
                    info[key] = value
        return pymet.SBMLReaction(unique_id=identifier, substrates=substrates,
                products=products, reversible=reaction.getReversible(),
                name=name, notes=info, **params)

