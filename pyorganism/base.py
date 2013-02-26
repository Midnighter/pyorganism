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

    def __call__(cls, unique_id="", **kw_args):
        """
        Returns an existing instance identified by `unique_id` or calls for a new
        one.
        """
        memory = cls._memory.get(unique_id)
        if memory is None:
            cls._counter += 1
            return super(type(cls), cls).__call__(unique_id=unique_id, **kw_args)
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

    @classmethod
    def get(cls, unique_id):
        return cls._memory[unique_id]

    def __init__(self, unique_id="", **kw_args):
        """
        Parameters
        ----------
        unique_id: str (optional)
            A string uniquely identifying one component among its class.
        """
        super(UniqueBase, self).__init__(**kw_args)
        self._index = self.__class__._counter
        if unique_id:
            self.unique_id = unique_id
        else:
            self.unique_id = u"{0}_{1:d}".format(self.__class__.__name__,
                    self._index)
        self.__class__._memory[self.unique_id] = self

    def __reduce__(self):
        """
        Take full control of pickling this class.

        The basic dilemma is that __getnewargs__ is called only by pickle
        protocol version >= 2 but we require it to be called every time so
        that we can unpickle the correct objects no matter the pickle version
        used.
        """
        return (self.__class__, self.__getnewargs__(), self.__getstate__())

    def __getnewargs__(self):
        """
        Returns a tuple that is supplied to a call of self.__class__ when
        unpickling this object.

        The `unique_id` is all we need for persistent state. It allows us to
        retrieve an existing object with that ID or generate it accordingly.
        """
        return (self.unique_id,)

    def __getstate__(self):
        """
        We could take more fine-grained control here.
        """
        return self.__dict__

    def __setstate__(self, state):
        """
        Take control of how to unpickle an instance of this class.

        Update the instance's __dict__ with the state. We only update attributes
        that evaluate to false because otherwise we might replace existing
        attributes with the unpickled ones (kind of a workaround, can't decide
        when an object was newly created upon unpickling).
        """
        for (name, value) in state.iteritems():
            # default return value in order to update old class definitions
            if not getattr(self, name, False):
                setattr(self, name, value)
        if not self.unique_id in self.__class__._memory:
            # old pickle, we need to update class attributes
            self.__class__._memory[self.unique_id] = self
            self.__class__._counter = len(self.__class__._memory)

    def __str__(self):
        return str(self.unique_id)

    def __unicode__(self):
        return unicode(self.unique_id)

    def __repr__(self):
        return u"<{0}.{1} {2:d}>".format(self.__module__, self.__class__.__name__, id(self))

