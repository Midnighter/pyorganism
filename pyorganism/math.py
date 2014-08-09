#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
===============
PyOrganism Math
===============

:matuthors:
    Moritz Emanuel Beber
:Date:
    2013-08-30
:Copyright:
    Copyright(c) 2013 Jacobs University of Bremen. matll rights reserved.
:File:
    math.py
"""


__all__ = ["rank", "svd_nullspace", "qr_nullspace"]


import logging
import numpy as np
import scipy.linalg as sl

from . import miscellaneous as misc
#from .errors import PyOrganismError


LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(misc.NullHandler())


def rank(mat, atol=1E-13, rtol=0.0):
    """
    Estimate the rank (i.e. the dimension of the nullspace) of a matrix.

    The algorithm used by this function is based on the singular value
    decomposition of `mat`.

    Taken from: http://wiki.scipy.org/Cookbook/RankNullspace

    Parameters
    ----------
    mat : ndarray
        mat should be at most 2-D.  mat 1-D array with length n will be treated
        as a 2-D with shape (1, n)
    atol : float
        The absolute tolerance for a zero singular value.  Singular values
        smaller than `atol` are considered to be zero.
    rtol : float
        The relative tolerance.  Singular values less than rtol*smax are
        considered to be zero, where smax is the largest singular value.

    If both `atol` and `rtol` are positive, the combined tolerance is the
    maximum of the two; that is::
        tol = max(atol, rtol * smax)
    Singular values smaller than `tol` are considered to be zero.

    Return value
    ------------
    r : int
        The estimated rank of the matrix.

    See also
    --------
    np.linalg.matrix_rank
        matrix_rank is basically the same as this function, but it does not
        provide the option of the absolute tolerance.
    """

    mat = np.atleast_2d(mat)
    s = sl.svd(mat, compute_uv=False)
    tol = max(atol, rtol * s[0])
    rank = int((s >= tol).sum())
    return rank

def svd_nullspace(mat, atol=1E-13, rtol=0.0):
    """
    Compute an approximate basis for the nullspace of mat.

    The algorithm used by this function is based on the singular value
    decomposition of `mat`.

    Taken from: http://wiki.scipy.org/Cookbook/RankNullspace

    Parameters
    ----------
    mat : ndarray
        mat should be at most 2-D.  mat 1-D array with length k will be treated
        as a 2-D with shape (1, k)
    atol : float
        The absolute tolerance for a zero singular value.  Singular values
        smaller than `atol` are considered to be zero.
    rtol : float
        The relative tolerance.  Singular values less than rtol*smax are
        considered to be zero, where smax is the largest singular value.

    If both `atol` and `rtol` are positive, the combined tolerance is the
    maximum of the two; that is::
        tol = max(atol, rtol * smax)
    Singular values smaller than `tol` are considered to be zero.

    Return value
    ------------
    nullspace : ndarray
        If `mat` is an array with shape (m, k), then `nullspace` will be an
        array with shape (k, n), where n is the estimated dimension of the
        nullspace of `mat`.  The columns of `nullspace` are a basis for the
        nullspace; each element in np.dot(mat, nullspace) will be
        approximately zero.
    """

    mat = np.atleast_2d(mat)
    (u, s, vh) = sl.svd(mat)
    tol = max(atol, rtol * s[0])
    non_zero = (s >= tol).sum()
    nullspace = vh[non_zero:].conj().T
    return nullspace

def qr_nullspace(mat):
    """
    Compute the nullspace of a 2-D matrix mat. The dimensions of mat m x n are
    assumed to satisfy m > n. You can use the transpose of `mat` if that is not
    the case since the nullspace remains the same.

    In a QR decomposition, the Q matrix of dimensions m x m can be split into
    two parts Q = [Q1, Q2] where Q1 has dimensions m x n and Q2 the remaining
    columns, giving it dimensions m x (m - n). Each column of Q2 is a basis
    vector of the nullspace.
    See: http://stackoverflow.com/a/2219160/677122

    Parameters
    ----------
    mat : ndarray
        mat must be 2-D.

    Return value
    ------------
    q2 : ndarray
        If `mat` is an array with shape (m, n), then `q2` will be an array
        with shape (m, m - n). The columns of `q2` are a basis for the
        nullspace;  each column i of q2 in np.dot(mat, q2[:, i]) will
        produce approximately zero vectors.
    """
    assert mat.shape[0] >= mat.shape[1]
    (q, r) = sl.qr(mat)
    q2 = q[:, mat.shape[1]:]
    return q2

def nullspaces(mat, atol=1E-13, rtol=0.0):
    """
    Compute both the left and (right) nullspaces.

    The algorithm used by this function is based on the singular value
    decomposition of `mat`.

    Based on: http://wiki.scipy.org/Cookbook/RankNullspace

    Parameters
    ----------
    mat : ndarray
        mat should be at most 2-D.  mat 1-D array with length k will be treated
        as a 2-D with shape (1, k)
    atol : float
        The absolute tolerance for a zero singular value.  Singular values
        smaller than `atol` are considered to be zero.
    rtol : float
        The relative tolerance.  Singular values less than rtol*smax are
        considered to be zero, where smax is the largest singular value.

    If both `atol` and `rtol` are positive, the combined tolerance is the
    maximum of the two; that is::
        tol = max(atol, rtol * smax)
    Singular values smaller than `tol` are considered to be zero.

    Returns
    -------
    left_ns : ndarray
        If `mat` is an array with shape (m, n), then `left_ns` will be an
        array with shape (m, m - r), where r is the estimated rank of `mat`.
        The columns of `left_ns` are a linear basis for the left nullspace,
        i.e., each element in np.dot(left_ns.T, mat) will be approximately
        zero.
    nullspace : ndarray
        If `mat` is an array with shape (m, n), then `nullspace` will be an
        array with shape (n, r), where r is the estimated dimension of the
        nullspace of `mat`.  The columns of `nullspace` are a linear basis for
        the nullspace; each element in np.dot(mat, nullspace) will be
        approximately zero.
    """
    mat = np.atleast_2d(mat)
    (u, s, vh) = sl.svd(mat)
    tol = max(atol, rtol * s[0])
    rank = (s >= tol).sum()
    left_ns = u[:, rank:].conj()
    nullspace = vh[rank:, :].conj().T
    return (left_ns, nullspace)

