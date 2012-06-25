#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
====================
Metabolic Components
====================

:Authors:
    Moritz Emanuel Beber
    Nikolaus Sonnenschein
:Date:
    2011-04-07
:Copyright:
    Copyright(c) 2011 Jacobs University of Bremen. All rights reserved.
:File:
    elements.py
"""


__all__ = ["BasicCompound", "BasicReaction", "BasicCompartment", "SBMLCompound",
        "SBMLCompartment", "SBMLCompartmentCompound", "SBMLReaction"]


import logging
import itertools

from .. import miscellaneous as misc
from ..errors import PyOrganismError
from ..base import UniqueBase


LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(misc.NullHandler())


class BasicCompound(UniqueBase):
    """
    The simplest form of representing a metabolic compound - just its class and
    name.
    """

    def __init__(self, name="", *args, **kw_args):
        """
        Parameters
        ----------
        name: str (optional)
            A string uniquely identifying the compound among its class.
        """
        super(BasicCompound, self).__init__(name, *args, **kw_args)


class BasicReaction(UniqueBase):
    """
    The simplest form of representing a biochemical reaction - just its class
    and name.

    Notes
    -----
    This is useful in metabolic networks where substrates and products are
    determined by the topology. Stoichiometric information is then stored on the
    links themselves. The only other information stored about the reaction is
    its reversibility.
    """

    def __init__(self, name="", reversible=False, *args, **kw_args):
        """
        Parameters
        ----------
        name: str (optional)
            A string uniquely identifying the reaction among its class.
        reversible: bool (optional)
            Reversibility information of the reaction.
        """
        super(BasicReaction, self).__init__(name, *args, **kw_args)
        self.reversible = reversible


class BasicCompartment(UniqueBase):
    """
    The simplest form of representing a cellular compartment - just its class
    and name.
    """

    def __init__(self, name="", *args, **kw_args):
        """
        Parameters
        ----------
        name: str
            A string uniquely identifying the compartment among its class.
        """
        super(BasicCompartment, self).__init__(name, *args, **kw_args)
        self._contained = set()

    def register(self, element):
        """
        Parameters
        ----------
        element: `BasicCompound`
            Compound that is found in this compartment.
        """
        self._contained.add(element)

    def __contains__(self, element):
        """
        Tests for the existance of `element` in this compartment.
        """
        if isinstance(element, BasicCompound):
            return element in self._contained
        elif isinstance(element, SBMLReaction):
            return all(cmpd in self._contained for cmpd in element)
        else:
            raise PyOrganismError(u"unrecognised metabolic component '{0}'",
                    element)


class BasicCompartmentCompound(BasicCompound):
    """
    A compartment specific compound.

    Often it is desirable to identify compounds on a per compartment basis, for
    example, in FBA experiments. This class is a simple container for both the
    compound instance that already exists and the compartment.
    """

    def __init__(self, name="", compound=None, compartment=None, *args, **kw_args):
        """
        Parameters
        ----------
        compound: BasicCompound
            An instance of BasicCompound that is then attached to a compartment.
        compartment: BasicCompartment
            An instance of BasicCompartment in which the compound is located.
        """
        super(BasicCompartmentCompound, self).__init__(name, *args,
                **kw_args)
        self.compound = compound
        self.compartment = compartment
        if not self.compartment is None:
            self.compartment.register(self.compound)

    def __getattr__(self, attr):
        return self.compound.__getattribute__(attr)


class SBMLCompartment(BasicCompartment):
    """
    A cellular compartment as defined per SBML standard.
    """

    def __init__(self, name="", outside=None, constant=True, suffix="",
            spatial_dimensions=None, size=None, units=None, *args, **kw_args):
        """
        Parameters
        ----------
        name: str
            A string uniquely identifying the compartment among its class.
        outside: str
            The name of the compartment that surrounds this one.
        constant: bool (optional)
            Determines whether the size attribute is allowed to change during
            model simulation.
        suffix: str (optional)
            A string appended to compounds for input/output.
        spatial_dimensions: int (optional)
            From 0 to 3, normal models have three dimensions.
        size: float (optional)
            The magnitude of the spatial_dimension in units.
        units: str (optional)
            A string identifying the unit in which size is measured.

        Notes
        -----
        The constant attribute is so far only kept for compatibility with SBML,
        it's not actually required. This behaviour may change in future.
        """
        super(SBMLCompartment, self).__init__(name, *args, **kw_args)
        self.outside = misc.convert(outside, SBMLCompartment)
        self.constant = bool(constant)
        self.suffix = suffix
        self.spatial_dimensions = misc.convert(spatial_dimensions, int)
        self.size = size
        self.units = units

    def __contains__(self, item):
        if isinstance(item, SBMLReaction):
            return all(cmpd.compartment == self for cmpd in item)
        elif isinstance(item, SBMLCompartmentCompound):
            return item.compartment == self
        else:
            raise PyOrganismError("unrecognised metabolic component '%s'", item)


class SBMLCompound(BasicCompound):
    """
    A molecular compound as defined per SBML standard.
    """

    def __init__(self, identifier, extended_name="", formula=None, in_chl=None,
            in_chl_key=None, smiles=None, charge=None, mass=None, *args, **kw_args):
        """
        Parameters
        ----------
        identifier: str
            A shorthand string uniquely identifying the compound among its class.
        extended_name: str
            A string uniquely identifying the compound.
        formula: str (optional)
            Molecular formula as a simple string, e.g., C6H12O6.
        in_chl: str (optional)
            An IUPAC compliant identifier in InChl format.
        in_chl_key: int (optional)
            A hashed key of the InChl string.
        smiles: str (optional)
            A SMILES representation of the compound.
        charge: int (optional)
            Electric charge on the compound (may be pH dependent).
        mass: float (optional)
            A unit-less magnitude determining the mass of the compound.
        """
        super(SBMLCompound, self).__init__(identifier, *args, **kw_args)
        self.identifier = identifier
        self.extended_name = extended_name
        self.formula = formula
        self.in_chl = in_chl
        self.in_chl_key = in_chl_key
        self.smiles = smiles
        self.charge = misc.convert(charge, int)
        self.mass = misc.convert(mass, float)

    def __contains__(self, element):
        """
        Checks for the existance of an atomic element in the compound.
        """
        raise NotImplementedError


class SBMLCompartmentCompound(BasicCompartmentCompound):
    """
    A compartment specific compound.

    Often it is desirable to identify compounds on a per compartment basis, for
    example, in FBA experiments. This class is a simple container for both the
    compound instance that already exists and the compartment.
    """

    def __init__(self, name="", compound=None, compartment=None, *args, **kw_args):
        """
        Parameters
        ----------
        compound: SBMLCompound
            An instance of SBMLCompound that is then attached to a compartment.
        compartment: SBMLCompartment
            An instance of SBMLCompartment in which the compound is located.
        """
        super(SBMLCompartmentCompound, self).__init__(name, compound,
                compartment, *args, **kw_args)

class SBMLReaction(BasicReaction):
    """
    A biochemical reaction as defined per SBML standard.
    """

    def __init__(self, identifier, substrates, products, reversible=False,
            extended_name="", synonyms=None, rate_constant=None,
            lower_bound=None, upper_bound=None, objective_coefficient=None,
            flux_value=None, reduced_cost=None, notes=False, *args, **kw_args):
        """
        Parameters
        ----------
        identifier: str
            A shorthand string uniquely identifying the reaction among its class.
        substrates: dict
            A map from the reaction educts to the aboslute value of their
            stoichiometric factors in the reaction.
        products: dict
            A map from the reaction products to the aboslute value of their
            stoichiometric factors in the reaction.
        reversible: bool (optional)
            Whether this reaction is known to occur in both directions in an
            organism.
        extended_name: str
            A string uniquely identifying the reaction.
        synonyms: str (optional)
            Additional identifiers of the reaction.
        rate_constant: float (optional)
            Unit-less specifier of the rate of the reaction at model conditions.
        notes: dict (optional)
            Additional notes, for example, from parsing an SBML model.
        """
        super(SBMLReaction, self).__init__(identifier, reversible, *args,
                **kw_args)
        self.identifier = identifier
        self.extended_name = extended_name
        self.substrates = substrates
        self.products = products
        self.reversible = bool(reversible)
        self.synonyms = synonyms
        self.rate_constant = misc.convert(rate_constant, float)
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound
        self.objective_coefficient = misc.convert(objective_coefficient, float)
        self.flux_value = flux_value
        self.notes = misc.convert(notes, dict, dict())
#        self._consistency_check()

    def __contains__(self, compound):
        """
        Parameters
        ----------
        compound: SBMLCompound
            A compound instance whose participation in the reaction is tested.
        """
        return (compound in self.substrates or compound in self.products)

    def __len__(self):
        return len(self.substrates) + len(self.products)

    def full_form(self):
        """
        Returns
        -------
        str:
            A string representation of the reaction, e.g., '2 A + 4 B -> 1 C'
            or '2 A + 4 B <=> 1 C' for a reversible reaction.
        """

        def util(compounds):
            for cmpd in compounds:
                yield str(abs(self.stoichiometric_coefficient(cmpd)))
                yield str(cmpd)
                if not (cmpd == compounds[-1]):
                    yield "+"
        rxn = ["%s:" % str(self.identifier)]
        rxn.extend([e for e in util(self.substrates.keys())])
        if self.reversible:
            rxn.append("<=>")
        else:
            rxn.append("->")
        rxn.extend([e for e in util(self.products.keys())])
        return " ".join(rxn)

    def compounds(self, coefficients=False):
        """
        Returns
        -------
        iterator:
            An iterator over all compounds partaking in the reaction.
        coefficients: bool (optional)
            Specifies whether the returned iterator should contain pairs of
            compounds with stoichiometric coefficients.
        """
        if coefficients:
            educts_iter = ((cmpd, -factor) for (cmpd, factor) in
                    self.substrates.iteritems())
            products_iter = ((cmpd, factor) for (cmpd, factor) in
                    self.products.iteritems())
            return itertools.chain(educts_iter, products_iter)
        else:
            return itertools.chain(self.substrates.iterkeys(),
                    self.products.iterkeys())

    def stoichiometric_coefficient(self, compound):
        """
        Parameters
        ----------
        compound: SBMLCompound
            A compound instance whose stoichiometric coefficient is sought for.

        Returns
        -------
        float:
            The stoichiometric coefficient of a compound in the reaction.
            Coefficients of substrates are negative.

        Exceptions
        ----------
        KeyError:
            In case compound is not part of the reaction.
        """
        if compound in self.substrates:
            return -self.substrates[compound]
        elif compound in self.products:
            return self.products[compound]
        else:
            raise KeyError("'{0}' is not participating in reaction"\
                    " '{1}'".format(str(compound), str(self)))

    def _consistency_check(self):
        """
        Asserts some basic consistency of the SBMLReaction instance.
        With enough meta data (SBMLCompound formula, charge, or mass)
        stoichiometric balancing is checked.

        Exceptions
        ----------
        AssertionError:
            In case any of the given conditions are not true.
        """
        # elemental balancing
        if all(cmpd.formula for cmpd in self.substrates.keys() +
                self.products.keys()):
            pass # not implemented yet
        # mass balancing
        if all(cmpd.mass for cmpd in self.substrates.keys() +
                self.products.keys()):
            assert sum(cmpd.mass * coeff for (cmpd, coeff) in
                    self.substrates.iteritems()) == sum(cmpd.mass * coeff
                    for (cmpd, coeff) in self.products.iteritems()),\
                    "There is a mass imbalance in reaction '{0}'".format(\
                    self.identifier)
        # charge balancing
        if all(cmpd.charge for cmpd in self.substrates.keys() +
                self.products.keys()):
            assert sum(cmpd.charge * coeff for (cmpd, coeff) in
                    self.substrates.iteritems()) == sum(cmpd.charge * coeff
                    for (cmpd, coeff) in self.products.iteritems()),\
                    "There is a charge imbalance in reaction '{0}'".format(\
                    self.identifier)

    def is_substrate(self, compound):
        """
        Parameters
        ----------
        compound: SBMLCompound
            A compound instance whose status as educt or product is queried.
        """
        return compound in self.substrates

