#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
==========
Exceptions
==========

:Author:
    Moritz Emanuel Beber
:Date:
    2012-05-22
:Copyright:
    Copyright(c) 2012 Jacobs University of Bremen. All rights reserved.
:File:
    errors.py
"""


__all__ = ["PyOrganismError"]


class PyOrganismError(StandardError):
    """
    An error for all exceptions that occur in the usage of the pymetabolism
    package.
    """

    def __init__(self, msg, *args, **kw_args):
        """
        Parameters
        ----------
        msg: str
            An unformatted string, i.e., it may contain multiple string format
            markers.

        Notes
        -----
        A variable number of arguments may be passed. They will all be used to
        format msg. So take care that the number and type of additional
        arguments matches the format markers in msg.

        Examples
        --------
        >>> err = PyMetabolismError("It's too {0} outside!", "rainy")
        >>> print(err)
        It's too rainy outside!
        >>> print(err.errno)
        1
        """
        super(PyOrganismError, self).__init__(msg, *args)
        self.errno = kw_args.get("errno", 1)
        if isinstance(msg, str):
            self.strerror = msg.format(*args, **kw_args)
        else:
            self.strerror = ""

    def __str__(self):
        return str(self.strerror)

    def __unicode__(self):
        return unicode(self.strerror)

