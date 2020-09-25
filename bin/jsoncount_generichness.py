#! /usr/bin/env python3
# -*- coding: ascii -*-

"""
Transform gene count matrices in derivatives.

Usage:
   countmatrixdb2generichness.py [-h] [-v] COUNTMATRIX_JSON GENERICHNESS_NREADS GENERICHNESS_ITER OUTPUT_PREFIX

Options:
  -h, --help
  -v, --version

"""

__author__ = "Matthieu Pichaud"
__version__ = "0.1"

from docopt import docopt
import sys
import os
import pandas
import numpy
import sqlite3
import json
import datetime
import gzip
import random

def unlist(ll):
    o = []
    for l in ll:
        o += l
    return o

def vector_generichness(v, downsizing_target=1000000, downsizing_iter=30):
    ## Create a list with as many item copies as counts for the item
    v_ = unlist([[k] * n for k, n in v])
    # print(len(set(v_)))
    if len(set(v_)) < downsizing_target:
        return -1
    ## Compute gene_richness
    generichness = []
    for i in range(downsizing_iter):
        generichness.append(len(set(random.sample(v_, downsizing_target))))
    return generichness

def countmatrix_generichness_onesample(
    countmatrix_df,
    outfile,
    downsizing_target=1000000,
    downsizing_iter=30,
    additionalinfo_d={},
):

    gr = vector_generichness(
        countmatrix_df.reset_index().values.tolist(), downsizing_target, downsizing_iter
    )

    data_d = {
        "generated_by": f"{os.path.abspath(sys.argv[0])} version:{__version__}",
        "date": str(datetime.datetime.now()),
        "data": [
            {
                "id": countmatrix_df.columns.tolist()[0],
                "downsizing_target": downsizing_target,
                "downsizing_iter": downsizing_iter,
                "mean": numpy.mean(gr),
                "std": numpy.std(gr),
            }
        ],
    }

    for k, v in additionalinfo_d.items():
        data_d[k] = v

    ## Export
    if outfile[-3:] != ".gz":
        outfile = outfile + ".gz"
    with gzip.open(outfile, "wt") as fo:
        json.dump(data_d, fo, indent=4)
    print(f"# - json file created: {outfile}")

    return outfile

if __name__ == "__main__":

    ## Collect arguments
    arguments = docopt(
        __doc__, version="{} version:{}".format(os.path.basename(__file__), __version__)
    )

    countmatrixjson = arguments["COUNTMATRIX_JSON"]
    generichness_nreads = int(arguments["GENERICHNESS_NREADS"])
    generichness_iter = int(arguments["GENERICHNESS_ITER"])
    output_prefix = arguments["OUTPUT_PREFIX"]

    ## Get started
    print(f"# {os.path.basename(__file__)}")
    print("# - load data")

    countmatrix_json = json.load(gzip.open(countmatrixjson, "rt"))
    countmatrix_df = pandas.json_normalize(countmatrix_json["data"])
    countmatrix_df = countmatrix_df.set_index(["id"])
    countmatrix_df = countmatrix_df.transpose()

    ## Collect parameters
    additionalinfo_d = {
        "reference": countmatrix_json["reference"],
        "reference_md5": countmatrix_json.get("reference_md5", "NA"),
    }

    del countmatrix_json

    ## Run
    countmatrix_generichness_onesample(
        countmatrix_df,
        f"{output_prefix}.generichness.json.gz",
        downsizing_target=generichness_nreads,
        downsizing_iter=generichness_iter,
        additionalinfo_d=additionalinfo_d,
    )
