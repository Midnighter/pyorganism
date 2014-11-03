# -*- coding: utf-8 -*-


"""
==================
PyOrganism Utility
==================

:Authors:
    Moritz Emanuel Beber
:Date:
    2014-11-01
:Copyright:
    Copyright(c) 2014 Jacobs University of Bremen. All rights reserved.
:File:
    utils.py
"""


__all__ = ["version_from_path"]


import logging
#import multiprocessing
#import signal
from os.path import (basename, dirname)
#from Queue import Empty

from . import miscellaneous as misc
#from .errors import PyOrganismError


LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(misc.NullHandler())


def version_from_path(path):
    ver = basename(path)
    if not ver:
        ver = basename(dirname(path))
    return ver

#def ignorant_worker(job_queue, result_queue):
#    signal.signal(signal.SIGINT, signal.SIG_IGN)
#    while True:
#        try:
#            job = job_queue.get(block=False)
#            result_queue.put(do_work())
#        except Empty:
#            pass
#
#
#class MultiPool(object):
#
#    def __init__(self, num_proc, **kw_args):
#        super(MultiPool, self).__init__(**kw_args)
#        self.num_proc = int(num_proc)
#        self.jobs = multiprocessing.Queue()
#        self.results = multiprocessing.Queue()
#        self.workers = [multiprocessing.Process(target=ignorant_worker,
#                args=(self.jobs, self.results)) for i in range(self.num_proc)]
#        for w in self.workers:
#            w.start()
