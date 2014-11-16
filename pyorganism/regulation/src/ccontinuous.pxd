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
#    ccontinuous.pxd
#
#.. |c| unicode:: U+A9
#"""


cdef extern from "continuous.h":
    double abs_control(int *sources, int *targets,
            const int num_links, double *expression)
    double difference_control(int *sources, int *targets,
            const int num_links, double *expression)
    double abs_difference_control(int *sources, int *targets,
            const int num_links, double *expression)
    double functional_control(int *sources, int *targets, int *function,
            const int num_links, double *expression)
    double functional_comparison(int *sources, int *targets, int *function,
            const int num_links, double *rate)
    double delayed_functional_control(int *sources, int *targets, int *function,
            const int num_links, double *expression, double *delta_expression)

