# -*- coding: utf-8 -*-


"""
===========================
PyOrganism Parallel Utility
===========================

:Authors:
    Moritz Emanuel Beber
:Date:
    2013-08-22
:Copyright:
    Copyright(c) 2013 Jacobs University of Bremen. All rights reserved.
:File:
    utils.py
"""


import time


def iter_async_results(async_results, timeout=1E-03):
    """
    This performs the same as iterating over a map_async result but is useful
    for a list of AsyncResult instances that were generated from calling
    view.apply many times. (This is useful for calling a remote function
    without arguments.)
    """
    pending = set(async_results)
    while pending:
        ready = set(ar for ar in pending if ar.ready())
        if not ready:
            time.sleep(timeout)
            continue
        # update pending to exclude those that just finished
        pending.difference_update(ready)
        for ar in ready:
            # we know these are done, so don't worry about blocking
            yield ar.get()

def multi_apply(view, function, num, *args, **kw_args):
    results = [view.apply(function, *args, block=False, ordered=False, **kw_args) for i in xrange(num)]
    return iter_async_results(results)

