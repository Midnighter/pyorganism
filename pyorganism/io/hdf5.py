# -*- coding: utf-8 -*-


"""
==========================
Input/Output of HDF5 Files
==========================

:Authors:
    Moritz Emanuel Beber
:Date:
    2013-12-10
:Copyright:
    Copyright(c) 2013 Jacobs University of Bremen. All rights reserved.
:File:
    hdf5.py
"""


__all__ = ["ResultManager"]


import os
import logging
import uuid

import numpy
import tables

from .. import miscellaneous as misc


LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(misc.NullHandler())

UUID_LENGTH = 35


class ExpressionData(tables.IsDescription):
    """
    """
    feature = tables.StringCol(6)
    value = tables.Float64Col() # normalized or not, see table description


class SimpleData(tables.IsDescription):
    """
    """
    value = tables.Float64Col() # either random sample or jackknife result


class ControlData(tables.IsDescription):
    """
    Other things to record (as attributes or in data structure):
        * link to TRN (version)
        * link to GPN (version)
        * link to metabolic network
        * sample tables named same as uuid
        * jackknife analysis named same as uuid
    """
    uuid = tables.StringCol(UUID_LENGTH) # however long uuids are
    db_version = tables.StringCol(4) # xx.x
    control_type = tables.StringCol(12) # analog or digital or metabolic
    continuous = tables.BoolCol() # discrete or continuous
    control = tables.Float64Col() # analog or digital control or metabolic coherence
    ctc = tables.Float64Col()
    robustness_mean = tables.Float64Col()
    robustness_std = tables.Float64Col()
    description = tables.StringCol(62)


class ResultManager(object):

    def __init__(self, filename, key="/", **kw_args):
        super(ResultManager, self).__init__(**kw_args)
        self.key = key
        if not os.path.exists(filename):
            self._setup(filename)
        else:
            self.h5_file = tables.open_file(filename, mode="a")
            self.root = self.h5_file._get_node(self.key)
            self.samples = self.root.samples
            self.robustness = self.root.robustness
            self.control = self.root.control

    def append(self, version, control_type, continuous, control_strength, ctc,
            description, samples=None, robustness=None):
        unique_id = "sim" + str(uuid.uuid4()).replace("-", "")
        unique_id = unique_id[:UUID_LENGTH]
        if samples is not None:
            table = self.h5_file.create_table(self.samples, unique_id, SimpleData,
                    expectedrows=len(samples))
            sample = table.row
            for item in samples:
                sample["value"] = item
                sample.append()
            table.flush()
        if robustness is not None:
            table = self.h5_file.create_table(self.robustness, unique_id, SimpleData,
                    expectedrows=len(robustness))
            sample = table.row
            for item in robustness:
                sample["value"] = item
                sample.append()
            table.flush()
        row = self.control.row
        row["uuid"] = unique_id
        row["db_version"] = version
        row["control_type"] = control_type
        row["continuous"] = continuous
        row["control"] = control_strength
        row["ctc"] = ctc
        row["description"] = description
        if robustness is not None:
            mask = numpy.isfinite(robustness)
            row["robustness_mean"] = robustness[mask].mean()
            row["robustness_std"] = robustness[mask].std(ddof=1)
        else:
            row["robustness_mean"] = numpy.nan
            row["robustness_std"] = numpy.nan
        row.append()
        self.control.flush()

    def _setup(self, filename):
        self.h5_file = tables.open_file(filename, mode="w",
                title="Control Analysis Data")
        components = self.key.split("/")
        self.h5_file.create_group("/", components[1])
        for i in range(2, len(components)):
            self.h5_file.create_group("/".join(components[:i]), components[i])
        self.root = self.h5_file._get_node(self.key)
        self.samples = self.h5_file.create_group(self.root, "samples", title="Null Model Samples")
        self.robustness = self.h5_file.create_group(self.root, "robustness",
                title="Robustness Analysis Results")
        self.control = self.h5_file.create_table(self.root, "control",
                ControlData, title="Summary table for all control analysis results.")

    def finalize(self):
        self.h5_file.close()

