# -*- coding: utf-8 -*-
#cython: boundscheck=False, wraparound=False


#"""
#===============
#Boolean Updates
#===============
#
#:Author:
#    Moritz Emanuel Beber
#:Date:
#    2014-11-01
#:Copyright:
#    Copyright |c| 2014, Jacobs University Bremen gGmbH, all rights reserved.
#:File:
#    _continuous.pyx
#
#.. |c| unicode:: U+A9
#"""


cimport ccontinuous as con


def abs_control_timeline(int[:] sources, int[:] targets, int num_links,
        double[:, ::1] expression, double[:] control, int num_points, int dim,
        **kw_args):
    cdef int i
    cdef double *expr_ptr = &expression[0, 0]
    for i in range(num_points):
        control[i] = con.abs_control(&sources[0], &targets[0], num_links,
                expr_ptr + (i * dim))

def difference_control_timeline(int[:] sources, int[:] targets, int num_links,
        double[:, ::1] expression, double[:] control, int num_points, int dim,
        **kw_args):
    cdef int i
    cdef double *expr_ptr = &expression[0, 0]
    for i in range(num_points):
        control[i] = con.difference_control(&sources[0], &targets[0], num_links,
                expr_ptr + (i * dim))

def abs_difference_control_timeline(int[:] sources, int[:] targets, int num_links,
        double[:, ::1] expression, double[:] control, int num_points, int dim,
        **kw_args):
    cdef int i
    cdef double *expr_ptr = &expression[0, 0]
    for i in range(num_points):
        control[i] = con.abs_difference_control(&sources[0], &targets[0],
                num_links, expr_ptr + (i * dim))

def functional_control_timeline(int[:] sources, int[:] targets,
        int[:] functions, int num_links, double[:, ::1] expression,
        double[:] control, int num_points, int dim, **kw_args):
    cdef int i
    cdef double *expr_ptr = &expression[0, 0]
    for i in range(num_points):
        control[i] = con.functional_control(&sources[0], &targets[0],
                &functions[0], num_links, expr_ptr + (i * dim))

def functional_comparison_timeline(int[:] sources, int[:] targets,
        int[:] functions, int num_links, double[:, ::1] expression, double[:] control,
        int num_points, int dim, **kw_args):
    cdef int i
    cdef double *expr_ptr = &expression[0, 0]
    for i in range(num_points):
        control[i] = con.functional_comparison(&sources[0], &targets[0],
                &functions[0], num_links, expr_ptr + (i * dim))

def delayed_functional_control_timeline(int[:] sources, int[:] targets,
        int[:] functions, int num_links, double[:, ::1] expression,
        double[:] control, int num_points, int dim, int delay, **kw_args):
    cdef int i
    cdef double *expr_ptr = &expression[0, 0]
    for i in range(num_points):
        control[i] = con.delayed_functional_control(&sources[0], &targets[0],
                &functions[0], num_links, expr_ptr + i * dim,
                expr_ptr + ((i + delay) % num_points) * dim)

