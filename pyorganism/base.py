#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
=======================
PyOrganism Base Classes
=======================

:Authors:
    Moritz Emanuel Beber
:Date:
    2012-06-03
:Copyright:
    Copyright(c) 2012 Jacobs University of Bremen. All rights reserved.
:File:
    base.py
"""


__all__ = ["UniqueBase"]


import logging

from . import miscellaneous as misc


LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(misc.NullHandler())


class MetaBase(type):
    """
    Metaclass for UniqueBase.

    This metaclass is one of two possible solutions to having a per-class
    dictionary that stores existing instances of that class. The per-class
    aspect is tricky because _memory is a mutable object. It is thus shared
    among UniqueBase and all its subclasses.

    Shared state of the _memory variable is avoided with this metaclass by
    creating each class with a different empty dict.

    Overriding the __call__ method ensures that existing instances (identified
    by the unique_id) are not re-initialised with different arguments.
    """
    def __new__(mcls, cls_name, cls_bases, cls_dct):
        """
        Adds a unique `dict` to each class.
        """
        cls_dct["_memory"] = dict()
        return super(MetaBase, mcls).__new__(mcls, cls_name, cls_bases, cls_dct)

    def __call__(cls, unique_id="", *args, **kw_args):
        """
        Returns an existing instance identified by `unique_id` or calls for a new
        one.
        """
        memory = cls._memory.get(unique_id)
        if memory is None:
            cls._counter += 1
            return super(type(cls), cls).__call__(unique_id, *args, **kw_args)
        else:
            return memory


class UniqueBase(object):
    """
    Base class for all objects that should be unique based on an identifier.

    Notes
    -----
    Each instance of this class or its subclasses is uniquely identified and
    stored by its identifier.  Instantiating the same class with the same
    identifier will simply yield the original instance.

    Since mutable class attributes are shared among the class and its subclasses
    there are two things to keep in mind when subclassing:
        1. A class specific storage of instances can be achieved by either using
        the provided metaclass or by giving each subclass its own attribute of
        `_memory = dict()`.
        2. The metaclass approach is more elegant since overriding the class'
        `__call__` method avoids re-initialising any existing instance with
        different arguments. It avoids having to put `if unique_id in _memory:
        return` into each subclasses' `__init__` method.

    The class attribute `_counter` is immutable and thus the state is not shared
    among subclasses but could be modified by the metaclass in the same way as
    `_memory` otherwise.

    Warning
    -------
    Subclasses of `UniqueBase` must not override `__call__` unless you know
    exactly what you're doing.

    Examples
    --------
    """

    # metaclass adds mutable subclass-specific attribute
    __metaclass__ = MetaBase
    # immutable class attribute is subclass-specific automatically
    _counter = 0

    def __init__(self, unique_id="", **kw_args):
        """
        Parameters
        ----------
        unique_id: str (optional)
            A string uniquely identifying one component among its class.
        """
        super(UniqueBase, self).__init__(**kw_args)
        # reading class attribute _counter
        self._index = self._counter
        if unique_id:
            self.unique_id = unique_id
        else:
            self.unique_id = u"{0}_{1:d}".format(self.__class__.__name__, self._index)
        # assigning to class attribute _memory
        self.__class__._memory[self.unique_id] = self

    def __str__(self):
        return str(self.unique_id)

    def __unicode__(self):
        return unicode(self.unique_id)

    def __repr__(self):
        return u"<{0}.{1} {2:d}>".format(self.__module__, self.__class__.__name__, id(self))

