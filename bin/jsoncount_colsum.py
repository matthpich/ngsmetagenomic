#! /usr/bin/env python3
# -*- coding: ascii -*-

"""
Colsum gene count matrix..

Usage:
   jsoncount_colsum.py [-h] [-v] [-n <RICHNESS_NREADS>] COUNTMATRIX_JSON OUTPUT_PREFIX

Options:
  -n <RICHNESS_NREADS>      Number of reads to downsize the matrix to before richness is computed [default: 0].
  -h, --help                Show this screen.
  -v, --version             Show version.

"""

__author__ = "Matthieu Pichaud"
__version__ = "0.1"

from docopt import docopt
import sys
import os
import pandas
import json
import gzip
import numpy
import sqlite3
import datetime
import hashlib
import shutil

def md5(fname):
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()

def countmatrix_colsum(countmatrix_df, outfile, richness_nreads=0, additionalinfo_d={}):

    print(f"# - countmatrix_norm", flush=True)
    colsum_df = countmatrix_df.sum(axis=0)
    colsum_d = colsum_df.to_dict()

    data_d = {
        "generated_by": f"{os.path.abspath(sys.argv[0])} version:{__version__}",
        "date": str(datetime.datetime.now()),
        "data": [],
    }

    if richness_nreads > 0:
        print(f"# - generichness like", flush=True)
        scaledcountmatrix_df = (
            ((countmatrix_df.div(colsum_df, axis=1) * richness_nreads) >= 0.5) * 1
        ).sum(axis=0)
        richness_d = scaledcountmatrix_df.to_dict()
        data_d["richness"] = richness_d
        data_d["richness_Nreads"] = richness_nreads

    for k in colsum_d.keys():
        data_d["data"].append({"id": k, "colsum": colsum_d.get(k, -1)})

    for k, v in additionalinfo_d.items():
        data_d[k] = v

    ## Export
    print("# - export", flush=True)
    if outfile[-3:] != ".gz":
        outfile = outfile + ".gz"
    with gzip.open(outfile, "wt") as fo:
        json.dump(data_d, fo, indent=4)
    print(f"# - json file created: {outfile}", flush=True)

    return outfile

if __name__ == "__main__":

    ## Collect arguments
    arguments = docopt(
        __doc__, version="{} version:{}".format(os.path.basename(__file__), __version__)
    )

    countmatrixjson = arguments["COUNTMATRIX_JSON"]
    output_prefix = arguments["OUTPUT_PREFIX"]
    richness_nreads = int(arguments["-n"])

    ## Get started
    print(f"# {os.path.basename(__file__)}")
    print("# - load data")

    countmatrix_json = json.load(gzip.open(countmatrixjson, "rt"))
    countmatrix_df = pandas.io.json.json_normalize(countmatrix_json["data"])
    countmatrix_df = countmatrix_df.set_index(["id"])
    countmatrix_df = countmatrix_df.transpose()

    ## Collect parameters
    additionalinfo_d = {
        "reference": countmatrix_json.get("reference", "NA"),
        "reference_md5": countmatrix_json.get("reference_md5", "NA"),
    }

    ## Run
    countmatrix_colsum(
        countmatrix_df,
        f"{output_prefix}.colsum.json.gz",
        richness_nreads=richness_nreads,
        additionalinfo_d=additionalinfo_d,
    )
