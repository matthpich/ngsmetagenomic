#! /usr/bin/env python3
# -*- coding: ascii -*-

"""
Transform MSP count matrix in presence/absence matrix

Usage:
   jsonMSP_modulepresenceabsence.py [-h] [-v] MSPCOUNTMATRIX_JSONGZ MSPDESCR_FN OUTPUT_PREFIX

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
import json
import gzip
import datetime
import hashlib
import shutil

def md5(fname):
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()

def mspcountmatrix_presenceabsence(
    mspcountmatrix_df, mspdescr_df, outfile, additionalinfo_d={}
):

    print(f"# - mspcountmatrix_presenceabsence", flush=True)

    mspmodule_all_presenceabsence_l = []

    for index, template_core in mspcountmatrix_df.iloc[
        mspcountmatrix_df.index.get_level_values("module_name") == "core", :
    ].iterrows():

        ##
        ##     Expected
        ##

        Nsamples = len(template_core)

        template_alpha = (
            mspdescr_df.iloc[
                mspdescr_df.index.get_level_values("msp_name") == index[0], :
            ]
            .reset_index()
            .drop("msp_name", axis=1)
            .set_index("module_name")
            .loc[:, "alpha"]
        )
        Naccessory = len(template_alpha)
        if Naccessory == 0:
            continue

        temp_core = pandas.DataFrame([template_core] * Naccessory)
        temp_alpha = pandas.DataFrame([template_alpha] * Nsamples)

        temp_alpha_columns = temp_alpha.columns
        temp_core_columns = temp_core.columns

        temp_core = temp_core.T
        temp_core.columns = temp_alpha_columns
        temp_core = temp_core.T
        temp_alpha.index = temp_core_columns
        temp_alpha = temp_alpha.T

        mspmodule_expected_df = numpy.log10(1 + temp_core) + temp_alpha
        mspmodule_expected_mask_df = ((numpy.log10(1 + temp_core) + temp_alpha) > 1) * 1

        ##
        ##     Observed
        ##

        mspmodule_observed_df = numpy.log10(
            1
            + mspcountmatrix_df.iloc[
                (mspcountmatrix_df.index.get_level_values("msp_name") == index[0])
                & (mspcountmatrix_df.index.get_level_values("module_name") != "core"),
                :,
            ].reset_index(level=0, drop=True)
        )

        ##
        ##     Score
        ##

        mspmodule_score = (
            mspmodule_observed_df
            / mspmodule_expected_df.loc[mspmodule_observed_df.index, :]
        )
        mspmodule_score = mspmodule_score[mspmodule_expected_mask_df == 1]

        mspmodule_all_presenceabsence_l.append(
            pandas.concat([mspmodule_score], keys=[index[0]], names=["msp_name"])
        )

    presabs_d = (
        pandas.concat(mspmodule_all_presenceabsence_l)
        .dropna(axis=0, how="all")
        .to_dict()
    )
    data_l = []
    for k, v in presabs_d.items():
        data_l.append(v)
        data_l[-1] = {
            f"{k[0]}:{k[1]}": _ for k, _ in data_l[-1].items() if not numpy.isnan(_)
        }
        data_l[-1]["id"] = k

    data_d = {
        "generated_by": f"{os.path.abspath(sys.argv[0])} version:{__version__}",
        "date": str(datetime.datetime.now()),
        "data": data_l,  # [derivativecount_d]
    }

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

    mspcountmatrixjsongz = arguments["MSPCOUNTMATRIX_JSONGZ"]
    mspdescr_fn = arguments["MSPDESCR_FN"]
    output_prefix = arguments["OUTPUT_PREFIX"]

    ## Get started
    print(f"# {os.path.basename(__file__)}")
    print("# - load data")

    ## Load MSP data
    mspcount_json_d = json.load(gzip.open(mspcountmatrixjsongz))
    mspcountmatrix_df = (
        pandas.DataFrame(mspcount_json_d["data"]).set_index("id").T.fillna(0)
    )
    mspcountmatrix_df = mspcountmatrix_df.reset_index()
    mspcountmatrix_df["msp_name"], mspcountmatrix_df["module_name"] = (
        mspcountmatrix_df["index"].str.split("|", 1).str
    )
    mspcountmatrix_df = mspcountmatrix_df.set_index(["msp_name", "module_name"]).drop(
        "index", axis=1
    )

    ## Load MSP descr data
    mspdescr_df = pandas.read_csv(mspdescr_fn, sep="\t", header=0, index_col=[0, 1])

    ## Collect parameters
    additionalinfo_d = {
        "reference": mspcount_json_d["reference"],
        "reference_md5": mspcount_json_d["reference_md5"],
        "derivative_db": mspcount_json_d["derivative_db"],
        "derivative_db_md5": mspcount_json_d["derivative_db_md5"],
        "mspdescription_fn": mspdescr_fn,
        "mspdescription_fn.md5": md5(mspdescr_fn),
    }

    del mspcount_json_d

    ## Run
    mspcountmatrix_presenceabsence(
        mspcountmatrix_df,
        mspdescr_df,
        f"{output_prefix}.MSP.modulepresenceabsence.json.gz",
        additionalinfo_d=additionalinfo_d,
    )
