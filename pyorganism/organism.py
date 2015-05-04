# -*- coding: utf-8 -*-


from __future__ import (absolute_import, unicode_literals)


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

from builtins import str
from future.utils import python_2_unicode_compatible

from . import miscellaneous as misc


LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(misc.NullHandler())

OPTIONS = misc.OptionsManager.get_instance()


@python_2_unicode_compatible
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

    def __init__(self, name, **kw_args):
        """
        Parameters
        ----------
        name: str
            The name of the organism. Should be unique but that's not a
            requirement.
        """
        super(Organism, self).__init__(**kw_args)
        self.name = str(name)
        self.genes = None
        self.trn = None
        self.gpn = None
        self.go = None
        self.couplons = None
        self.metabolism = None
        self.metabolic_network = None
        self.activity = dict()
        self.significant = dict()

    def __str__(self):
        return self.name

