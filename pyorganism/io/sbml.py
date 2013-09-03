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


__all__ = ["SBMLParser", "ARABIDOPSIS_THALIANA", "ESCHERICHIA_COLI_textbook",
        "ESCHERICHIA_COLI_iAF1260", "ESCHERICHIA_COLI_iJO1366",
        "SACCHAROMICES_CEREVISIAE_iND750", "HOMO_SAPIENS_Recon1"]


import os
import logging
#import re

from .. import miscellaneous as misc
from ..metabolism import elements as pyel
from ..errors import PyOrganismError

libsbml = misc.load_module("libsbml", "SBML", "http://sbml.org/Software/libSBML")


LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(misc.NullHandler())

OPTIONS = misc.OptionsManager.get_instance()

ARABIDOPSIS_THALIANA =  {"Peroxisome": "_x", "Cytoplasm": "_c",
        "Plastid": "_p", "Golgi_appratus": "_g",
        "Vacuole": "_v", "Mitochondrion": "_m",
        "Endoplasmic_reticulum": "_r"}

ESCHERICHIA_COLI_textbook = {"C_c": "_c", "C_e": "_e"}
ESCHERICHIA_COLI_iAF1260 = {"C_c": "_c", "C_p": "_p", "C_e": "_e"}
ESCHERICHIA_COLI_iJO1366 = {"c": "_c", "p": "_p", "e": "_e"}

HOMO_SAPIENS_Recon1 = {"C_x": "_x", "C_e": "_e", "C_g": "_g", "C_c": "_c",
        "C_n": "_n", "C_r": "_r", "C_m": "_m", "C_l": "_l"}

SACCHAROMICES_CEREVISIAE_iND750 = {"C_x": "_x", "C_e": "_e", "C_g": "_g",
        "C_c": "_c", "C_v": "_v", "C_n": "_n", "C_r": "_r", "C_m": "_m"}

class SBMLParser(object):
    """
    A class implementing methods for parsing a SBML model
    """

    def __init__(self, **kw_args):
        super(SBMLParser, self).__init__(**kw_args)
        self.document = None
        self.compartment_ids = None
        self.compound_ids = None
        self.reaction_ids = None
#        self.reaction_property = re.compile("([a-z]+[a-z_0-9]*)$", re.UNICODE)
        self._key = None
        self._value = None

    def __call__(self, model, **kw_args):
        """
        Parse a document in SBML format.
        """
        if isinstance(model, basestring):
            if os.path.exists(model):
                # might want to read file on our own to determine encoding
                self.document = libsbml.readSBMLFromFile(model)
            else:
                self.document = libsbml.readSBMLFromString(model)
        elif isinstance(model, file):
            self.document = libsbml.readSBMLFromString(model.read())
        num_err = self.document.getNumErrors()
        fatal = False
        if num_err > 0:
            for i in xrange(num_err):
                sbml_err = self.document.getError(i)
                LOGGER.warn(sbml_err.getShortMessage())
                if sbml_err.getSeverity() == libsbml.LIBSBML_SEV_FATAL:
                    fatal = True
        if fatal:
            raise PyOrganismError("fatal error in parsing SBML document")
        return self.parse()

    def parse(self):
        from ..metabolism.systems import MetabolicSystem
        self.compartment_ids = dict()
        self.compound_ids = dict()
        self.reaction_ids = dict()
        self._model = self.document.getModel()
        for compartment in self._model.getListOfCompartments():
            self._parse_compartment(compartment)
        LOGGER.info("%d compartments", len(self.compartment_ids))
        for cmpd in self._model.getListOfSpecies():
            self._parse_species(cmpd)
        LOGGER.info("%d (compartmentalized) compounds", len(self.compound_ids))
        for rxn in self._model.getListOfReactions():
            self._parse_reaction(rxn)
        LOGGER.info("%d reactions", len(self.reaction_ids))
        return MetabolicSystem(compartments=self.compartment_ids.values(),
                compounds=self.compound_ids.values(),
                reactions=self.reaction_ids.values())

    def _parse_notes_paragraph(self, par, info, note_sep=":"):
        if not (self._key is None or self._value is None):
            info[self._key] = self._value
            self._key = None
            self._value = None
        content = par.getChild(0).toString().strip()
        if content.endswith(note_sep):
            if not self._key is None:
                info[self._key] = None
            self._key = content[:-1].strip().lower().replace(" ", "_")
            return
        tmp = content.split(note_sep, 1)
        if len(tmp) == 1:
            if not self._value is None:
                LOGGER.warn("replacing note value: '%s'", self._value)
            self._value = content
        elif len(tmp) == 2:
            if not self._key is None:
                info[self._key] = None
            elif not self._value is None:
                LOGGER.warn("replacing note value: '%s'", self._value)
            self._key = tmp[0].strip().lower().replace(" ", "_")
            self._value = tmp[1].strip()
        else:
            LOGGER.warn("unrecognized notes format: '%s'",
                    par.getChild(0).toString())

    def _parse_notes(self, element):
        # other information contained in notes
        self._key = None
        self._value = None
        info = dict()
        notes = element.getNotes()
        if notes is None:
            return info
        for i in range(notes.getNumChildren()):
            # some notes sections contain just <p> tags, others are enclosed in
            # additional <body>...</body> section
            node = notes.getChild(i)
            tag = node.toString()
            if tag == "<p>":
                self._parse_notes_paragraph(node, info)
            elif tag == "<body>":
                for j in range(node.getNumChildren()):
                    child = node.getChild(j)
                    child_tag = child.toString()
                    if child_tag == "<p>":
                        self._parse_notes_paragraph(child, info)
                    else:
                        LOGGER.warn("unexpected child XML tag '%s'", child_tag)
            elif tag == "<html>":
                for j in range(node.getNumChildren()):
                    child = node.getChild(j)
                    child_tag = child.toString()
                    if child_tag == "<p>":
                        self._parse_notes_paragraph(child, info)
                    else:
                        LOGGER.warn("unexpected child XML tag '%s'", child_tag)
            else:
                LOGGER.warn("unexpected XML tag '%s'", tag)
        return info

    def _parse_compartment(self, compartment):
        unique = compartment.getId()
        suffix = OPTIONS.compartment_suffixes[unique]
        self.compartment_ids[compartment.getId()] =  pyel.SBMLCompartment(
                unique_id=unique,
                name=compartment.getName(),
                outside=compartment.getOutside(),
                constant=compartment.getConstant(),
                spatial_dimensions=compartment.getSpatialDimensions(),
                size=compartment.getSize(), units=compartment.getUnits(),
                suffix=suffix)

    def _parse_species(self, compound):
        """
        Able to parse entries from getListOfSpecies

        @todo: Check for meta information and parse if available
        """
        compartment = self.compartment_ids.get(compound.getCompartment())
        identifier = compound.getId()
        identifier = identifier.replace("_DASH_", "-")
#        identifier = identifier.replace("_LPAREN_", "(")
#        identifier = identifier.replace("_RPAREN_", ")")
        if identifier.startswith(OPTIONS.compound_prefix):
            identifier = identifier[len(OPTIONS.compound_prefix):]
        if compartment is None:
            for (unique, suffix) in OPTIONS.compartment_suffixes.iteritems():
                if identifier.endswith(suffix):
                    identifier = identifier[:-len(suffix)]
                    compartment = self.compartment_ids[unique]
                    break
        else:
            identifier = identifier[:-len(compartment.suffix)]
        name = compound.getName()
        info = self._parse_notes(compound)
        cmpd = pyel.SBMLCompound(unique_id=identifier, name=name,
                charge=compound.getCharge(), notes=info)
        if compartment is None:
            self.compound_ids[compound.getId()] = cmpd
        else:
            self.compound_ids[compound.getId()] = pyel.SBMLCompartmentCompound(
                    unique_id=identifier + compartment.suffix,
                    compound=cmpd, compartment=compartment)

    def _strip_reaction_id(self, name):
        identifier = name
        identifier = name.replace("_DASH_", "-")
#        identifier = identifier.replace("_LPAREN_", "(")
#        identifier = identifier.replace("_RPAREN_", ")")
        if identifier.startswith(OPTIONS.reaction_prefix):
            identifier = identifier[len(OPTIONS.reaction_prefix):]
        if identifier.endswith(OPTIONS.reversible_suffix):
            identifier = identifier[:-len(OPTIONS.reversible_suffix)]
        return identifier

    def _parse_reaction(self, reaction):
        """Able to parse entries from getListOfReactions"""
        identifier = self._strip_reaction_id(reaction.getId())
#        if OPTIONS.exchange_reaction in identifier:
#            # TODO
#            for cmpd in
#                (self._parse_species(model.getSpecies(elem.getSpecies())) for
#                        elem in reaction.getListOfReactions()):
        name = reaction.getName()
        # parse additional reaction parameters
        try:
            parameters = reaction.getKineticLaw().getListOfParameters()
            params = dict((param.getId().lower(), param.getValue())\
                    for param in parameters)
        except AttributeError:
            params = dict()
            LOGGER.debug("reaction '%s' has no kinetic parameters",
                    reaction.getId())
        # substrates' stoichiometry
        substrates = dict((self.compound_ids[elem.getSpecies()],
                abs(elem.getStoichiometry()))\
                for elem in reaction.getListOfReactants())
        # products' stoichiometry
        products = dict((self.compound_ids[elem.getSpecies()],
                abs(elem.getStoichiometry()))\
                for elem in reaction.getListOfProducts())
        # other information contained in notes
        info = self._parse_notes(reaction)
        # reaction properties and distinguishing suffixes for various
        # compartments cannot be separated easily but suffixes are necessary
#        mobj = self.reaction_property.search(identifier)
#        if mobj:
#            info["properties"] = mobj.group(1)
#            identifier = identifier[:mobj.start(1)]
        self.reaction_ids[reaction.getId()] = pyel.SBMLReaction(
                unique_id=identifier,
                substrates=substrates,
                products=products,
                reversible=reaction.getReversible(),
                name=name, notes=info, **params)

