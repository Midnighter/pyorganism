#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
==========
PyOrganism
==========

:Authors:
    Moritz Emanuel Beber
:Date:
    2012-05-22
:Copyright:
    Copyright(c) 2012 Jacobs University of Bremen. All rights reserved.
:File:
    organism.py
"""


__all__ = ["Organism"]


import logging

from .. import miscellaneous as misc
from ..errors import PyOrganismError


LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(misc.NullHandler())

OPTIONS = misc.OptionsManager.get_instance()


class Organism(object):
    """
    A representation of a living organism with multiple layers of organisation.

    As many layers of organisation as are available or desired may be included
    in the `Organism` object.

    Notes
    -----

    Examples
    --------
    """

    def __init__(self, name):
        """
        Parameters
        ----------
        name: str
            The name of the organism. Should be unique but that's not a
            requirement.
        """
        object.__init__(self)
        self.trn = None
        self.gpn = None
        self.go = None
        self.couplons = None
        self.metabolism = None

    def digital_control(self, active):
        """
        Compute the digital control from the effective TRN.

        Parameters
        ----------
        active: iterable
            An iterable with names of genes that are active in a specific
            condition.

        Warning
        -------
        Unknown gene names are silently ignored.
        """
        assert self.trn
        subnet = self.trn.subgraph(active)
        controlled = sum(1 for (node, deg) in subnet.degree_iter() if deg > 0)
        return controlled / float(subnet.order())

    def analog_control(self, active):
        """
        Compute the digital control from the effective GPN.

        Parameters
        ----------
        active: iterable
            An iterable with names of genes that are active in a specific
            condition.

        Warning
        -------
        Unknown gene names are silently ignored.
        """
        assert self.gpn
        subnet = self.gpn.subgraph(active)
        controlled = sum(1 for (node, deg) in subnet.degree_iter() if deg > 0)
        return controlled / float(subnet.order())

