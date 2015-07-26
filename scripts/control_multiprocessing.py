#!/usr/bin/env python
# -*- coding: utf-8 -*-


from __future__ import (absolute_import, unicode_literals, division,
        print_function)


import sys
import os
import logging
import argparse
import multiprocessing
from logging.config import dictConfig

import numpy as np
from future.utils import raise_
from sqlalchemy import create_engine
from sqlalchemy.orm import joinedload
from progressbar import (ProgressBar, Timer, SimpleProgress, Bar, Percentage, ETA)

import pyorganism as pyorg
import pyorganism.regulation as pyreg
import pyorganism.io.models as pymodels


LOGGER = logging.getLogger()
LOGGER.addHandler(logging.StreamHandler())


def check_path(path):
    if not os.path.isfile(path):
        raise IOError("file does not exist '{path}'".format(path=path))

def string_pair(key, value, level, sep=": "):
    return "{0}{1}{2}".format(str(key), sep, str(value)).rjust(level)

def print_dict(dic, level=0):
    message = list()
    for (key, value) in dic.iteritems():
        if isinstance(value, dict):
            message.append(string_pair(key, "{", level))
            message.extend(print_dict(value, level + 2))
            message.append(string_pair("", "}", level, sep=""))
        elif isinstance(value, list):
            message.append(string_pair(key, "[...]", level))
        elif isinstance(value, set):
            message.append(string_pair(key, "{...}", level))
        else:
            message.append(string_pair(key, value, level))
    return message


##############################################################################
# Discrete
##############################################################################


def load_discrete(organism, config):
    LOGGER.info("Loading differentially expressed genes:")
    experiments = config["experiments"]
    for (filename, name, reader, extra) in zip(experiments["paths"],
            experiments["names"], experiments["readers"],
            experiments["args"]):
        reader_func = getattr(pymicro, reader)
        path = os.path.join(experiments["base"], filename)
        check_path(path)
        LOGGER.info("  %s: '%s'", name, path)
        organism.activity[name] = reader_func(path)

def simple_discrete(control_type, df, name2gene):
    eligible = df["name"].notnull()
    active = [name2gene[name] for name in df["name"][eligible]\
            if name2gene[name] is not None]
    LOGGER.info("  mapped %d/%d active genes (%.2f%%)", len(active), len(df),
            len(active) / len(df) * 100.0)
    results = dict()
    results["gene"] = active
    if control_type == "digital":
        results["gene"] = pyreg.active_genes_and_tf(active)
    results["tu"] = pyreg.active_tu(active)
    results["operon"] = pyreg.active_operons(active)
    return results

def ratio_discrete(control_type, df, name2gene):
    eligible = df["name"].notnull()
    up = (df["ratio (A/B)"] > 1.0) & eligible
    down = (df["ratio (A/B)"] < 1.0) & eligible
    up_reg = [name2gene[name] for name in df["name"][up]\
            if name2gene[name] is not None]
    total = up.sum()
    LOGGER.info("  mapped %d/%d up-regulated genes (%.2f%%)", len(up_reg), total,
            len(up_reg) / total * 100.0)
    down_reg = [name2gene[name] for name in df["name"][down]\
            if name2gene[name] is not None]
    total = down.sum()
    LOGGER.info("  mapped %d/%d down-regulated genes (%.2f%%)", len(down_reg),
            total, len(down_reg) / total * 100.0)
    results = dict()
    results["gene"] = dict()
    results["gene"]["up"] = up_reg
    results["gene"]["down"] = down_reg
    results["tu"] = dict()
    results["tu"]["up"] = pyreg.active_tu(up_reg)
    results["tu"]["down"] = pyreg.active_tu(down_reg)
    results["operon"] = dict()
    results["operon"]["up"] = pyreg.active_operons(up_reg)
    results["operon"]["down"] = pyreg.active_operons(down_reg)
    return results

def discrete_jobs(organism, config, *args):
    LOGGER.info("Generating discrete job specifications:")
    jobs = list()
    analysis = config["analysis"]
    for version in config["versions"]:
        for (cntrl_name, experiments, setups,
                control, ctc, measures, random_num, robustness_num,
                rob_extra, projections) in zip(analysis["control_types"],
                analysis["experimental_sets"], analysis["experimental_setups"],
                analysis["control"], analysis["ctc"],
                analysis["measures"], analysis["random_num"],
                analysis["robustness_num"], analysis["robustness_args"],
                config["network"]["projections"]):
            for method in ctc:
                for basis in projections:
                    for ms_name in measures:
                        for (exp_name, exp_setup) in zip(experiments, setups):
                            if exp_setup == "ratio_discrete":
                                for direction in ["up", "down"]:
                                    spec = dict()
                                    spec["version"] = version
                                    spec["continuous"] = config["continuous"]
                                    spec["control_type"] = cntrl_name
                                    spec["experiment"] = exp_name
                                    spec["projection"] = basis
                                    spec["setup"] = exp_setup
                                    spec["direction"] = direction
                                    spec["control"] = control
                                    spec["ctc"] = method
                                    spec["measure"] = ms_name
                                    spec["random_num"] = random_num
                                    spec["robustness_num"] = robustness_num
                                    spec["robustness_args"] = rob_extra
                                    jobs.append(spec)
                            else:
                                spec = dict()
                                spec["version"] = version
                                spec["continuous"] = config["continuous"]
                                spec["control_type"] = cntrl_name
                                spec["experiment"] = exp_name
                                spec["projection"] = basis
                                spec["setup"] = exp_setup
                                spec["control"] = control
                                spec["ctc"] = method
                                spec["measure"] = ms_name
                                spec["random_num"] = random_num
                                spec["robustness_num"] = robustness_num
                                spec["robustness_args"] = rob_extra
                                jobs.append(spec)
    LOGGER.info("  %d jobs", len(jobs))
    return jobs

def discrete_worker(spec):
    LOGGER.debug(spec)
    version = spec["version"]
    cntrl_type = spec["control_type"]
    global_vars = globals()
    control = getattr(pyreg, spec["control"])
    ctc = getattr(pyreg, spec["ctc"])
    measure = getattr(pyreg, spec["measure"])
    net = global_vars["networks"][version][cntrl_type][spec["projection"]]
    if cntrl_type == "analog":
        active = global_vars["prepared"][version][cntrl_type][spec["experiment"]][spec["projection"]][spec["direction"]]
    else:
        active = global_vars["prepared"][version][cntrl_type][spec["experiment"]][spec["projection"]]
    LOGGER.debug(len(active))
    effective = pyreg.effective_network(net, active)
    res_cntrl = control(effective, measure=measure)
    (res_ctc, samples) = ctc(effective, net, measure=measure,
            random_num=spec["random_num"], return_sample=True)
    return (spec, res_cntrl, res_ctc, samples)

def discrete_result(manager, spec, res_cntrl, res_ctc, samples):
    manager.append(version=spec["version"], control_type=spec["control_type"],
            continuous=spec["continuous"], strain=spec["experiment"],
            projection=spec["projection"], setup=spec["setup"], control_strength=res_cntrl,
            control_method=spec["control"], ctc=res_ctc, ctc_method=spec["ctc"],
            measure=spec["measure"], direction=spec.get("direction", None),
            samples=samples)

def main_discrete(args):
    bar = ProgressBar(maxval=len(tasks), widgets=[Timer(), " ",
            SimpleProgress(), " ", Percentage(), " ", Bar(), " ",
            ETA()]).start()
    pass

##############################################################################
# Continuous
##############################################################################


def normed(session, experiment):
    series = pymodels.Expression.load_frame(session, experiment)
    return series.apply(pyorg.norm_zero2unity, axis=1, raw=True)

def shuffle_feature(session, experiment):
    series = normed(session, experiment)
    # shuffles rows relative to index (the features)
    np.random.shuffle(series.values)
    return series

def shuffle_series(session, experiment):
    series = normed(session, experiment)
    # shuffles columns relative to column names (time points)
    np.random.shuffle(series.values.T)
    return series

def shuffle_all(session, experiment):
    series = normed(session, experiment)
    # reshuffles all values (flat iterator over all values in the 2D array)
    np.random.shuffle(series.values.flat)
    return series

def continuous_exec(args):
    (control, measure, random_num, delay, job_id) = args
    if "comparison" in measure:
        points = control.expression.columns[1:]
    else:
        points = control.expression.columns
    # include points somehow in the results
    (z_scores, ctrl_scores, samples) = control.series_ctc(measure, random_num,
            delay)
    return (job_id, z_scores, ctrl_scores, samples, points)

def main_continuous(args):
    glbls = globals()
    engine = create_engine(args.engine)
    pymodels.Base.metadata.bind = engine
    pymodels.Session.configure(bind=engine)
    session = pymodels.Session()
    tasks = session.query(pymodels.Job).\
            options(joinedload("analysis"), joinedload("control"),
            joinedload("experiment")).filter(~pymodels.Job.complete).all()
    if len(tasks) == 0:
        LOGGER.warn("Nothing to do")
        return
    analysis_configs = {job.analysis for job in tasks}
    control_configs = {job.control for job in tasks}
    experiments = {job.experiment for job in tasks}
    preparations = {job.preparation for job in tasks}
    sampling = {job.sampling for job in tasks}
    projections = {job.projection for job in tasks}
    LOGGER.debug("%d analysis configurations", len(analysis_configs))
    LOGGER.debug("%d control configurations", len(control_configs))
    LOGGER.debug("%d experiments", len(experiments))
    LOGGER.debug("%d setup cases", len(preparations))
    LOGGER.debug("%d sampling methods", len(sampling))
    LOGGER.debug("%d network projections", len(projections))
    num_prep = len(analysis_configs) * len(control_configs) * len(experiments)\
            * len(preparations) * len(sampling) * len(projections)
    LOGGER.debug("%d total configurations", num_prep)
    LOGGER.info("Preparing Data")
    task_args = dict()
    bar = ProgressBar(maxval=num_prep, widgets=[Timer(), " ",
            SimpleProgress(), " ", Percentage(), " ", Bar(), " ",
            ETA()]).start()
    for anal in analysis_configs:
        LOGGER.debug(" %s:", anal.version)
        feature2node = pyorg.read_pickle(os.path.join(anal.objects, anal.map))
        for cntrl in control_configs:
            LOGGER.debug("  %s", cntrl.type)
            net = pyorg.read_pickle(os.path.join(anal.objects, cntrl.network))
            tu_net = pyreg.to_transcription_unit_based(net)
            op_net = pyreg.to_operon_based(net)
            for exp in experiments:
                LOGGER.debug("   %s", exp.strain)
                for prep in preparations:
                    LOGGER.debug("    %s", prep)
                    series = glbls[prep](session, exp)
                    for sampl in sampling:
                        LOGGER.debug("     %s", sampl)
                        for prj in projections:
                            LOGGER.debug("      %s", prj)
                            control = pyreg.ContinuousControl()
                            if prj == "tu":
                                control.setup(tu_net, series, feature2node, sampl)
                            elif prj == "operon":
                                control.setup(op_net, series, feature2node, sampl)
                            else:
                                control.setup(net, series, feature2node, sampl)
                            if cntrl.type == "analog":
                                control.from_gpn()
                            elif cntrl.type == "digital":
                                control.from_trn()
                            else:
                                raise ValueError("'{}'".format(cntrl.type))
                            task_args[(anal.id, cntrl.id, exp.id, prep, sampl, prj)] = control
                            bar += 1
    bar.finish()
    LOGGER.info("Running Jobs")
    tasks = [(task_args[(job.analysis.id, job.control.id, job.experiment.id,
        job.preparation, job.sampling, job.projection)],) + (job.measure,
        job.random_num, job.delay, job.id) for job in tasks]
    pool = multiprocessing.Pool(args.nproc)
    result_it = pool.imap_unordered(continuous_exec, tasks)
    bar = ProgressBar(maxval=len(tasks), widgets=[Timer(), " ",
            SimpleProgress(), " ", Percentage(), " ", Bar(), " ",
            ETA()]).start()
    for (job_id, z_scores, cntrl_scores, samples, points) in result_it:
        results = list()
        try:
            job = session.query(pymodels.Job).filter_by(id=job_id).one()
            for (i, name) in enumerate(points):
                res = pymodels.Result(control=cntrl_scores[i], ctc=z_scores[i],
                        point=name, job=job)
                session.add(res)
                results.append(res)
            job.complete = True
            session.commit()
        except Exception:
            session.rollback()
            bar += 1
            continue
        if job.selection > 0:
            try:
                for (i, res) in enumerate(results):
                    # use a more low-level insert for speed
                    session.execute(pymodels.RandomSample.__table__.insert(),
                            [{"control": val, "result_id": res.id}\
                            for val in np.random.choice(samples[i], job.selection,
                            replace=False)])
                session.commit()
            except Exception:
                session.rollback()
                bar += 1
                continue
        bar += 1
    bar.finish()
    session.close()


##############################################################################
# Main
##############################################################################


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=None)
    parser.add_argument("-v", "--version", action="version", version="0.1")
    parser.add_argument("--log-level", dest="log_level", default="INFO",
            help="Log level, i.e., DEBUG, INFO, WARN, ERROR, CRITICAL (default: %(default)s)")
    parser.add_argument("--encoding", dest="encoding", default="utf-8",
            help="File encoding to assume (default: %(default)s)")
    parser.add_argument("-n", "--nproc", dest="nproc",
            default=multiprocessing.cpu_count(), type=int,
            help="Number of processors to use (default: %(default)s)")
    parser.add_argument("engine",
            help="Database connection string, e.g., 'sqlite+pysqlite:///file.db'")
    subparsers = parser.add_subparsers(help="sub-command help")
    # discrete
    parser_discrete = subparsers.add_parser("discrete",
            help="Perform a discrete control analysis")
    parser_discrete.set_defaults(func=main_discrete)
    # continuous
    parser_continuous = subparsers.add_parser("continuous",
            help="Perform a continuous control analysis")
    parser_continuous.set_defaults(func=main_continuous)
    args = parser.parse_args()
    dictConfig({"version": 1, "incremental": True,
        "root": {"level": args.log_level}})
    if args.log_level != "DEBUG":
        logging.getLogger("pyorganism").setLevel(logging.WARN)
    try:
        sys.exit(args.func(args))
    except: # we want to catch everything
        (err, msg, trace) = sys.exc_info()
        # do something
        raise_(err, msg, trace)
    finally:
        logging.shutdown()

