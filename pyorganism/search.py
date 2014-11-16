# -*- coding: utf-8 -*-


"""
============================
Searching Object Collections
============================

:Authors:
    Moritz Emanuel Beber
:Date:
    2013-12-12
:Copyright:
    Copyright(c) 2013 Jacobs University of Bremen. All rights reserved.
:File:
    search.py
"""


__all__ = ["FindObject"]


import logging
import warnings

import numpy as np

from . import miscellaneous as misc


LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(misc.NullHandler())

try:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        from fuzzywuzzy import process
        from fuzzywuzzy import fuzz
except ImportError:
    LOGGER.warn("fuzzy search requires 'fuzzywuzzy'"\
            " (https://github.com/seatgeek/fuzzywuzzy)")


class FindObject(object):
    def __init__(self, collection, attr=None, default="", targets=None,
            indeces=None, **kw_args):
        """
        Sets up an instance that makes one attribute of a collection of objects
        searcheable or allows searching of an arbitrary target vector where
        indeces provided link back to the collection.

        The object can be initialised either with an attribute name (plus
        optionally a default value if the attribute doesn't exist) or with a
        pre-computed vector `targets` to be searched. The attributes can have
        any value that is comparable. In none string cases, a different default
        value should be passed.
        If a searchable vector is passed, each item should correspond to the
        position of an object in `collection`, otherwise another index vector
        must be passed that links the two.
        """
        super(FindObject, self).__init__(**kw_args)
        self.collection = list(collection)
        if attr is not None:
            self.targets = np.array([getattr(item, attr, default)\
                    for item in self.collection])
            self.indeces = self.targets.argsort()
            self.targets = self.targets[self.indeces]
        elif targets is not None:
            self.targets = np.array(targets)
            if indeces is None:
                self.indeces = np.arange(len(collection), dtype=int)
            else:
                self.indeces = np.array(indeces, dtype=int)
            assert len(self.targets) == len(self.indeces), "search targets is"\
                    " longer than indeces, please provide"
            tmp_i = self.targets.argsort()
            self.targets = self.targets[tmp_i]
            self.indeces = self.indeces[tmp_i]

    def __call__(self, value):
        return self.binary_search(value)

    def binary_search(self, value):
        index = np.searchsorted(self.targets, value)
        if index == len(self.targets) or self.targets[index] != value:
            raise IndexError("not found '%s'" % value)
        return self.collection[self.indeces[index]]

    def fuzzy_search(self, value, threshold=80, scorer=fuzz.QRatio):
        match = process.extractOne(value, self.targets, scorer=scorer, score_cutoff=threshold)
        if match is None:
            raise IndexError("'%s' not found" % value)
        LOGGER.debug("'%s' matches '%s' (%d%%)", value, match[0],
                match[1])
        return (self.binary_search(match[0]), match[0], match[1])

    def match_iter(self, value, operator=np.equal):
        return (self.collection[i] for i in self.indeces[operator(self.targets, value)])

