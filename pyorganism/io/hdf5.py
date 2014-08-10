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


import logging
import uuid

import numpy as np
import tables

from .. import miscellaneous as misc


LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(misc.NullHandler())

UUID_LENGTH = 32 # stripped dashes


class ExpressionData(tables.IsDescription):
    """
    """
    feature = tables.StringCol(6)
    value = tables.Float64Col() # normalized or not, see table description


class SampleData(tables.IsDescription):
    """
    Random sample data belonging to a particular analysis.
    """
    value = tables.Float64Col()
    study = tables.StringCol(UUID_LENGTH)


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
    strain = tables.StringCol(30) # experimental strain
    projection = tables.StringCol(30) # gene, TU, operon
    setup = tables.StringCol(30) # data setup
    description = tables.StringCol(30) # misc
    direction = tables.StringCol(4) # up or down (based on fold-change)
    control = tables.Float64Col() # analog or digital control or metabolic coherence
    control_method = tables.StringCol(60)
    ctc = tables.Float64Col()
    ctc_method = tables.StringCol(60)
    measure = tables.StringCol(60)
    robustness_mean = tables.Float64Col()
    robustness_std = tables.Float64Col()
    time = tables.UInt32Col()
    delay = tables.UInt32Col()


class ResultManager(object):

    def __init__(self, filename, key="/", **kw_args):
        super(ResultManager, self).__init__(**kw_args)
        self.key = key
        self.h5_file = tables.open_file(filename, mode="a",
                title="Control Analysis Data")
        self.samples = None
        self.robustness = None
        self.control = None
        self._setup()

    def append(self, version, control_type, continuous, strain, projection,
            setup, control_strength, control_method, ctc, ctc_method, measure,
            description=None, samples=None, robustness=None, direction=None,
            time=None, delay=None):
        unique_id = str(uuid.uuid4()).replace("-", "")
        unique_id = unique_id[:UUID_LENGTH]
        if samples is not None:
            row = self.samples.row
            for item in samples:
                row["value"] = item
                row["study"] = unique_id
                row.append()
            self.samples.flush()
        if robustness is not None:
            row = self.robustness.row
            for item in samples:
                row["value"] = item
                row["study"] = unique_id
                row.append()
            self.robustness.flush()
        row = self.control.row
        row["uuid"] = unique_id
        row["db_version"] = version
        row["control_type"] = control_type
        row["continuous"] = continuous
        row["strain"] = strain
        row["projection"] = projection
        row["setup"] = setup
        if description is not None:
            row["description"] = description
        if direction is not None:
            row["direction"] = direction
        row["control"] = control_strength
        row["control_method"] = control_method
        row["ctc"] = ctc
        row["ctc_method"] = ctc_method
        row["measure"] = measure
        if time is not None:
            row["time"] = int(time)
        if delay is not None:
            row["delay"] = int(delay)
        if robustness is not None:
            mask = np.isfinite(robustness)
            row["robustness_mean"] = robustness[mask].mean()
            row["robustness_std"] = robustness[mask].std(ddof=1)
        else:
            row["robustness_mean"] = np.nan
            row["robustness_std"] = np.nan
        row.append()
        self.control.flush()

    def _setup(self):
        root = self.h5_file.root
        components = self.key.split("/")
        for node in components:
            try:
                root = self.h5_file.get_node(root, node)
                if not isinstance(root, tables.Group):
                    raise IOError("given key '{key}' contains none-Group nodes".format(
                            key=self.key))
            except tables.NoSuchNodeError:
                root = self.h5_file.create_group(root, node)
        try:
            self.samples = root.samples
        except tables.NoSuchNodeError:
            self.samples = self.h5_file.create_table(root, "samples",
                    SampleData, title="Null Model Samples")
        try:
            self.robustness = root.robustness
        except tables.NoSuchNodeError:
            self.robustness = self.h5_file.create_table(root, "robustness",
                    SampleData, title="Robustness Analysis Results")
        try:
            self.control = root.control
        except tables.NoSuchNodeError:
            self.control = self.h5_file.create_table(root, "control",
                    ControlData, title="Summary for all control analysis results.")

    def finalize(self):
        self.h5_file.close()

