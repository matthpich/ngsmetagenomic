#! /usr/bin/env python3
# -*- coding: ascii -*-

"""
Transform gene count matrices in derivatives.

Usage:
   countmatrixdb2derivatives.py [-h] [-v] COUNTMATRIX_JSON DERIVATIVES_FN OUTPUT_PREFIX

Options:
  -h, --help
  -v, --version

"""

__author__ = "Matthieu Pichaud"
__version__ = "0.4"

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
import math

def md5(fname):
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()

def parse_derivativefile(derivatives_fn):
    """
    Parse the file that describes the multiple derivatives to compute.
    Format:
        line 1:     header
        line 2+:    derivative_name \t db_filepath \t db_table \t derivative_column1,derivative_column2,...
    """

    derivative_d = {}

    with open(derivatives_fn) as fi:
        header = 1
        for line in fi:
            if header:
                header = 0
                continue
            items = line.rstrip("\r\n").split("\t")
            derivative_d[items[0]] = [items[1], items[2], items[3].split(",")]

    return derivative_d

def countmatrix_derivate(
    countmatrix_df,
    derivative_db,
    derivative_table,
    derivative_cols,
    outfile,
    aggregate_fun=sum,
    additionalinfo_d={},
):

    print(f"# - countmatrix_derivate - {derivative_db}", flush=True)

    ## Get derivative data
    print("#   |_ get derivative data", flush=True)
    tmp_db = "tmp." + ".".join(derivative_db.split("/")[-1].split(".")[1:])
    shutil.copyfile(derivative_db, tmp_db)

    print(f"derivative_db: {derivative_db}")

    if derivative_db.endswith(".db"):

        print("# db-type")

        con = sqlite3.connect(tmp_db)
        cursor = con.cursor()  # create cursor
        query_select = f"SELECT * FROM  {derivative_table}"
        derivative_df = pandas.read_sql(query_select, con)
        con.close()
        os.remove(tmp_db)

        derivative_df = derivative_df.set_index(["gene_name"], drop=True)

    elif derivative_db.endswith((".txt", ".tsv", ".txt.gz", ".tsv.gz")):

        print("# txt- or tsv-type")

        derivative_df = pandas.read_csv(
            tmp_db, sep="\t", index_col="gene_name", header=0
        )
        os.remove(tmp_db)

    elif derivative_db.endswith((".csv", ".csv.gz")):

        print("# csv-type")

        derivative_df = pandas.read_csv(
            tmp_db, sep=",", index_col="gene_name", header=0
        )
        os.remove(tmp_db)

    ## ============================================================
    ##          main derivative
    ## ============================================================

    ## Join tables on gene_name
    print("#   |_ join tables on gene_name", flush=True)
    print("-", flush=True)
    print(derivative_df.shape, flush=True)
    print(countmatrix_df.shape, flush=True)
    derivativecount_df = derivative_df.join(countmatrix_df, how="inner")
    print("    done.", flush=True)
    ## Clean unused columns
    print("#   |_ clean unused column", flush=True)
    for cn in ["gene_id", "gene_name", "top30"]:
        while cn in list(derivativecount_df):
            derivativecount_df.drop([cn], axis=1, inplace=True)
    print("control - derivativecount_df shape", derivativecount_df.shape)

    ## Aggregate
    print("#   |_ aggregate counts", flush=True)
    derivativecount_df = (
        derivativecount_df.groupby(derivative_cols, as_index=True)
        .aggregate(sum)
        .replace(0, numpy.nan)
        .dropna(axis=0, how="all")
    )
    if isinstance(derivativecount_df.index, pandas.core.indexes.multi.MultiIndex):
        derivativecount_df.index = derivativecount_df.index.map("|".join)
    derivativecount_d = derivativecount_df.to_dict()
    data_l = []
    for k, v in derivativecount_d.items():
        vclean = {_[0]: _[1] for _ in v.items() if math.isnan(_[1]) == False}
        data_l.append(vclean)
        data_l[-1]["id"] = k
    # derivativecount_d["id"] = derivativecount_df.columns.tolist()[0]

    data_d = {
        "generated_by": f"{os.path.abspath(sys.argv[0])} version:{__version__}",
        "date": str(datetime.datetime.now()),
        "derivative_db": derivative_db,
        "derivative_db_md5": md5(derivative_db),
        "data": data_l,  # [derivativecount_d]
    }

    for k, v in additionalinfo_d.items():
        data_d[k] = v

    ## Export
    print("#   |_ export counts", flush=True)
    if outfile[-3:] != ".gz":
        outfile = outfile + ".gz"
    with gzip.open(outfile, "wt") as fo:
        json.dump(data_d, fo, indent=4)
    print(f"#   |_ json file created: {outfile}", flush=True)

    ## ============================================================
    ##          fingerprint
    ## ============================================================

    ## Join tables on gene_name
    print("#   |_ join tables on gene_name (fingerprint)", flush=True)
    derivativecount_df = derivative_df.join(
        countmatrix_df.mask(countmatrix_df > 0, 1), how="inner"
    )

    ## Clean unused columns
    print("#   |_ clean unused column (fingerprint)", flush=True)
    for cn in ["gene_id", "gene_name", "top30"]:
        while cn in list(derivativecount_df):
            derivativecount_df.drop([cn], axis=1, inplace=True)
    print("control - derivativecount_df shape", derivativecount_df.shape)

    ## Aggregate fingerprint
    print("#   |_ aggregate (fingerprint)", flush=True)
    derivativecount_df = (
        derivativecount_df.groupby(derivative_cols, as_index=True)
        .aggregate(sum)
        .replace(0, numpy.nan)
        .dropna(axis=0, how="all")
    )
    if isinstance(derivativecount_df.index, pandas.core.indexes.multi.MultiIndex):
        derivativecount_df.index = derivativecount_df.index.map("|".join)
    derivativecount_d = derivativecount_df.to_dict()
    data_l = []
    for k, v in derivativecount_d.items():
        vclean = {_[0]: _[1] for _ in v.items() if math.isnan(_[1]) == False}
        data_l.append(vclean)
        data_l[-1]["id"] = k

    data_d = {
        "generated_by": f"{os.path.abspath(sys.argv[0])} version:{__version__}",
        "date": str(datetime.datetime.now()),
        "derivative_db": derivative_db,
        "derivative_db_md5": md5(derivative_db),
        "data": data_l,  # [derivativecount_d]
    }

    for k, v in additionalinfo_d.items():
        data_d[k] = v

    ## Export
    print("#   |_ export (fingerprint)", flush=True)
    if outfile[-3:] != ".gz":
        outfile = outfile + ".fingerprint.gz"
    else:
        outfile = outfile.replace(".gz", ".fingerprint.gz")
    with gzip.open(outfile, "wt") as fo:
        json.dump(data_d, fo, indent=4)
    print(f"#   |_ json file created (fingerprint): {outfile}", flush=True)

    return outfile

if __name__ == "__main__":

    ## Collect arguments
    arguments = docopt(
        __doc__, version="{} version:{}".format(os.path.basename(__file__), __version__)
    )

    countmatrixjson = arguments["COUNTMATRIX_JSON"]
    derivatives_p = arguments["DERIVATIVES_FN"]
    output_prefix = arguments["OUTPUT_PREFIX"]

    ## Get started
    print(f"# {os.path.basename(__file__)}")
    print(f"# - load data {countmatrixjson}")

    countmatrix_json = json.load(gzip.open(countmatrixjson, "rt"))
    countmatrix_df = pandas.json_normalize(countmatrix_json["data"])
    countmatrix_df = countmatrix_df.set_index(["id"])
    countmatrix_df = countmatrix_df.transpose()

    ## Collect parameters
    additionalinfo_d = {
        "reference": countmatrix_json.get("reference", "NA"),
        "reference_md5": countmatrix_json.get("reference_md5", "NA"),
    }

    del countmatrix_json

    derivatives_d = parse_derivativefile(derivatives_p)
    for derivative_name in derivatives_d:
        print(f"\n# {derivative_name}")
        derivative_filechek = 0
        if os.path.exists(derivatives_d[derivative_name][0]):
            derivative_fn = derivatives_d[derivative_name][0]
            print(f"{derivative_fn} found.")
            derivative_filecheck = 1
        elif os.path.exists(os.path.basename(derivatives_d[derivative_name][0])):
            derivative_fn = os.path.basename(derivatives_d[derivative_name][0])
            print(f"{derivative_fn} found.")
            derivative_filecheck = 1
        else:
            print(f"{derivative_fn} not found!!!")
            derivative_filecheck = 0

        if derivative_filecheck:
            countmatrix_derivate(
                countmatrix_df,
                derivative_fn,
                derivatives_d[derivative_name][1],
                derivatives_d[derivative_name][2],
                f"{output_prefix}.{derivative_name}.json.gz",
                additionalinfo_d=additionalinfo_d,
            )
        else:
            print(f"{derivative_name} not found.")
