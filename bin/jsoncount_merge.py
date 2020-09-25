#!/usr/bin/env python3

"""
Merge count matrices.

Usage:
    jsoncount_merge.py [-h] [-v] <jsoncountA.gz> ...

Options:
  -h, --help
  -v, --version
"""

__author__ = "Matthieu Pichaud"
__version__ = "1.1"

from docopt import docopt
import os
import sys
import datetime
import gzip
import json

controlkeys = ["reference", "reference_md5"]

def json_deserialize(fn):
    with gzip.open(fn, "rt") as read_file:
        data_d = json.load(read_file)
    return data_d

def join_jsoncountmatrices(fn_l, output_path):

    # - collect data
    for ifn, fn in enumerate(fn_l):
        data_d = json_deserialize(fn)
        if ifn == 0:
            all_d = data_d
            all_d["generated_by"] = f"{sys.argv[0]} version:{__version__}"
            all_d["date"] = str(datetime.datetime.now())
        else:
            for controlkey in controlkeys:
                if data_d[controlkey] != all_d[controlkey]:
                    print("# Trying to merge incompatible json > abort")
                    quit()
            all_d["data"] = all_d["data"] + data_d["data"]

    # - write output
    with gzip.open(output_path, "wt") as fo:
        json.dump(all_d, fo, indent=4)

    return output_path

##
##      Main
##

if __name__ == "__main__":

    arguments = docopt(
        __doc__, version="{} version:{}".format(os.path.basename(__file__), __version__)
    )

    count_paths = arguments["<jsoncountA.gz>"]
    join_jsoncountmatrices(count_paths, "jsoncountmerge.tsv.gz")
