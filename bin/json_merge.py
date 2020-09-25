#!/usr/bin/env python3

"""
Merge json files

Usage:
    json_merge.py [-h] [-v] <json> ...

Options:
  -h, --help
  -v, --version
"""

__author__ = "Matthieu Pichaud"
__version__ = "1.0"

from docopt import docopt
import os
import sys
import datetime
import json


def json_deserialize(fn):
    with open(fn, "rt") as read_file:
        data_d = json.load(read_file)
    return data_d

def merge_dict(dict1, dict2):
    dictmerged = {}
    k1 = set(dict1.keys())
    k2 = set(dict2.keys())
    for k in k1 - k2: dictmerged[k] = dict1[k]
    for k in k2 - k1: dictmerged[k] = dict2[k]
    for k in k1.intersect(k2):
        if isinstance(dict1[k], dict):
		if isinstance(dict2[k], dict): dictmerged[k] = merge_dict(dict1[k], dict2[k])
		else: return error incompatible content
	
	
    



def join_json(fn_l, output_path):

    # - collect data
    for ifn, fn in enumerate(fn_l):
        data_d = json_deserialize(fn)
        if ifn == 0:
            all_d = data_d
        else:


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

    json_paths = arguments["<json>"]
    join_json(json_paths, "merged.json")
