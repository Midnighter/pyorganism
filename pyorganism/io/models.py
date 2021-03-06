# -*- coding: utf-8 -*-


"""
===============
Database Models
===============

:Authors:
    Moritz Emanuel Beber
:Date:
    2014-11-10
:Copyright:
    Copyright(c) 2014 Jacobs University of Bremen. All rights reserved.
:File:
    models.py
"""


__all__ = [
    "Experiment",
    "Expression",
    "AnalysisConfiguration",
    "KnockOut",
    "Job",
    "RandomSample",
    "Result"
]


import logging

from sqlalchemy import (Table, Column, ForeignKey, Integer, String, Sequence,
        Float, Boolean, select, and_)
from sqlalchemy.orm import (sessionmaker, relationship)
from sqlalchemy.ext.declarative import declarative_base
from pandas import DataFrame
from numpy import (nan, argsort)

from .. import miscellaneous as misc


LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(misc.NullHandler())


Session = sessionmaker()
Base = declarative_base()


class KnockOut(Base):
    __tablename__ = "knockout"
    id = Column(Integer, Sequence("knockout_id_seq"), primary_key=True)
    feature = Column(String(15)) # usually up to four letters


knockout_association = Table("ko_ass", Base.metadata,
    Column("knockout_id", Integer, ForeignKey("knockout.id")),
    Column("experiment_id", Integer, ForeignKey("experiment.id"))
)


class Expression(Base):
    __tablename__ = "expression"
    id = Column(Integer, Sequence("expression_id_seq"), primary_key=True)
    level = Column(Float)
    feature = Column(String(15)) # usually up to four letters
    point = Column(String(30)) # describe the point in a series
    experiment_id = Column(Integer, ForeignKey("experiment.id"))

    @classmethod
    def load_frame(cls, session, experiment):
        """
        Load part of the table into a well-formatted pandas.DataFrame.

        session can be any object with the execute method.
        """
        table = cls.__table__
        stmt = select([table.c.feature, table.c.point, table.c.level]).where(
                table.c.experiment_id == experiment.id)
        query = session.execute(stmt)
        df = DataFrame(iter(query), columns=query.keys())
        df.set_index(["feature", "point"], inplace=True)
        series = df.unstack()
        series.columns = series.columns.droplevel()
        # time points can become unsorted in database, sort them
        series = series.reindex_axis(
                series.columns[ argsort(series.columns.astype(int).values)],
                axis=1, copy=False)
        if experiment.knockouts is not None:
            series.loc[[ko.feature for ko in experiment.knockouts]] =  nan
        return series


class ExpressionRatio(Base):
    __tablename__ = "expressionratio"
    id = Column(Integer, Sequence("expressionratio_id_seq"), primary_key=True)
    control_id = Column(Integer, ForeignKey("expression.id"))
    control = relationship("Expression", uselist=False, foreign_keys=[control_id])
    treatment_id = Column(Integer, ForeignKey("expression.id"))
    treatment = relationship("Expression", uselist=False,
            foreign_keys=[treatment_id])
    log2_ratio = Column(Float)
    p_value = Column(Float)


class Experiment(Base):
    __tablename__ = "experiment"
    id = Column(Integer, Sequence("experiment_id_seq"), primary_key=True)
    expression = relationship("Expression", backref="experiment")
    knockouts = relationship("KnockOut", secondary=knockout_association)
    strain = Column(String(30))
    organism = Column(String())
    description = Column(String())


class AnalysisConfiguration(Base):
    __tablename__ = "analconfig"
    id = Column(Integer, Sequence("analconfig_id_seq"), primary_key=True)
    type = Column(String(10)) # continuous, discrete
    version = Column(String(5)) # usually just x.x
    objects = Column(String()) # db objects directory
    map = Column(String()) # e.g., feature2gene.pkl

class ControlConfiguration(Base):
    __tablename__ = "controlconfig"
    id = Column(Integer, Sequence("controlconfig_id_seq"), primary_key=True)
    type = Column(String(10)) # analog, digital, metabolic
    network = Column(String()) # e.g., trn.pkl
    window = Column(Integer) # only used for analog, usually 5000
    direction = Column(String(1)) # only used for analog, "+" or "-"

class Job(Base):
    __tablename__ = "job"
    id = Column(Integer, Sequence("job_id_seq"), primary_key=True)
    analysis_id = Column(Integer, ForeignKey("analconfig.id"))
    analysis = relationship("AnalysisConfiguration", uselist=False)
    control_id = Column(Integer, ForeignKey("controlconfig.id"))
    control = relationship("ControlConfiguration", uselist=False)
    experiment_id = Column(Integer, ForeignKey("experiment.id"))
    experiment = relationship("Experiment", uselist=False)
    results = relationship("Result", backref="job")
    preparation = Column(String()) # e.g., simple
    projection = Column(String(30)) # usually just gene, tu, operon
    sampling = Column(String(30)) # timeline, fork or fork-strand
    measure = Column(String()) # absolute, functional, fixed_tf and so on
    random_num = Column(Integer) # often 1E05
    delay = Column(Integer) # usually between 1 and 6 with delay measures
    selection = Column(Integer) # how many elements to store of the random_num
    complete = Column(Boolean, default=False, nullable=False)


class RandomSample(Base):
    __tablename__ = "sample"
    id = Column(Integer, Sequence("sample_id_seq"), primary_key=True)
    control = Column(Float)
    result_id = Column(Integer, ForeignKey("result.id"))

    @classmethod
    def load_frame(cls, session):
        """
        Load part of the table into a well-formatted pandas.DataFrame.

        session can be any object with the execute method.
        """
        sample = cls.__table__
        job = Job.__table__
        result = Result.__table__
        analysis = AnalysisConfiguration.__table__
        control = ControlConfiguration.__table__
        experiment = Experiment.__table__
        stmt = select([sample.c.id, sample.c.control,
                result.c.point, control.c.type, control.c.direction,
                experiment.c.strain, job.c.preparation, job.c.sampling,
                job.c.projection, job.c.measure, job.c.delay,
                analysis.c.version]).where(and_(
                sample.c.result_id == result.c.id,
                result.c.job_id == job.c.id,
                job.c.analysis_id == analysis.c.id,
                job.c.control_id == control.c.id,
                job.c.experiment_id == experiment.c.id))
        query = session.execute(stmt)
        df = DataFrame(iter(query), columns=query.keys())
        df.set_index("id", inplace=True)
        return df


class Result(Base):
    __tablename__ = "result"
    id = Column(Integer, Sequence("result_id_seq"), primary_key=True)
    job_id = Column(Integer, ForeignKey("job.id"))
    samples = relationship("RandomSample", backref="result")
    control = Column(Float)
    ctc = Column(Float)
    point = Column(String(30)) # describe the point in a series

    @classmethod
    def load_frame(cls, session):
        """
        Load part of the table into a well-formatted pandas.DataFrame.

        session can be any object with the execute method.
        """
        result = cls.__table__
        job = Job.__table__
        analysis = AnalysisConfiguration.__table__
        control = ControlConfiguration.__table__
        experiment = Experiment.__table__
        stmt = select([result.c.id, result.c.control, result.c.ctc,
                result.c.point, control.c.type, control.c.direction,
                experiment.c.strain, job.c.preparation, job.c.sampling,
                job.c.projection, job.c.measure, job.c.delay,
                analysis.c.version]).where(and_(
                result.c.job_id == job.c.id,
                job.c.analysis_id == analysis.c.id,
                job.c.control_id == control.c.id,
                job.c.experiment_id == experiment.c.id))
        query = session.execute(stmt)
        df = DataFrame(iter(query), columns=query.keys())
        df.set_index("id", inplace=True)
        return df

