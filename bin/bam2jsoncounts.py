#!/usr/bin/env python3
# -*- coding: ascii -*-

"""
Create a count matrix from a sorted_bam file.

Sort reads by name and count:
- if the 2 paired reads map the same target --> +1
- if the 2 paired reads map different targets --> +1,+1
- if one of the paired reads map a target --> +1

Usage:
    bam_map2count.py [--help] [--version] [--threads THREADS] SAMPLEID REFFNA BAM ...

Options:
  --help
  --version
  -t THREADS, --threads THREADS   number of threads used [default: 1]
"""

__author__ = "Matthieu Pichaud"
__version__ = "1.2"

from docopt import docopt
import subprocess
import os
import sys
import string
from collections import Counter
import json
import csv
import gzip
import hashlib
import datetime

##
##      Supporting functions
##

def execute(cmd):
    popen = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, universal_newlines=True, shell=True
    )
    for stdout_line in iter(popen.stdout.readline, ""):
        yield stdout_line
    popen.stdout.close()
    return_code = popen.wait()
    if return_code:
        raise subprocess.CalledProcessError(return_code, cmd)

def md5(fname):
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()

##
##      Main
##

if __name__ == "__main__":

    arguments = docopt(
        __doc__, version="{} version:{}".format(os.path.basename(__file__), __version__)
    )

    sampleid = arguments["SAMPLEID"]
    reffna = arguments["REFFNA"]
    bam = arguments["BAM"]
    Nthreads = arguments["--threads"]

    ##  Get counts
    count_d = {}
    neighbor_l = []
    for bam_ in bam:

        os.system(f"samtools sort -n -@{Nthreads} {bam_} > {bam_}.sorted.bam")
        cmd = f"samtools view -@{Nthreads} {bam_}.sorted.bam | grep -v 'XA:Z' | grep -v 'SA:Z'"

        readid_prev = ""
        hit_prev = ""

        for iline, line in enumerate(execute(cmd)):
            items = line.rstrip("\r\n").split("\t")
            if len(items) >= 2:
                readid = items[0]
                hit = items[2]
                if readid_prev == readid:
                    # increment the read1 hit
                    count_d[hit] = count_d.get(hit, 0) + 1
                    if hit_prev != hit:
                        # increment the read2 hit
                        count_d[hit_prev] = count_d.get(hit_prev, 0) + 1
                        # report neighbor
                        neighbor_l.append((min(hit_prev, hit), max(hit_prev, hit)))
                    # reset
                    readid_prev = ""
                    hit_prev = ""
                elif readid_prev != "":
                    # read is different from previous read, increment the previous read hit
                    count_d[hit_prev] = count_d.get(hit_prev, 0) + 1
                    readid_prev = readid
                    hit_prev = hit
                else:
                    # new set of read - store
                    readid_prev = readid
                    hit_prev = hit
        os.remove(f"{bam_}.sorted.bam")
    print("bam file parsed.")

    ##  Report pairs of sequences mapped by reads from the same pair
    c = Counter(neighbor_l)
    with gzip.open("bam2counts.neighbor.tsv.gz", "wt") as fo:
        fo_ = csv.writer(fo, delimiter="\t")
        for k, v in c.items():
            if v >= 5:
                fo_.writerow([k[0], k[1], v])
    print("neighbor information saved.")

    ##  Parse reffna
    with open(reffna, "rt") as fi:
        seq_d = {}
        seqcontent = ""
        for line in fi:
            line = line.rstrip("\r\n")
            if line[0] == ">":
                if len(seqcontent):
                    seq_d[seqid] = seqcontent
                seqid = line[1:]
                seqcontent = ""
            else:
                seqcontent += line
        if len(seqcontent):
            seq_d[seqid] = seqcontent

    ## Report counts
    data_d = {
        "generated_by": f"{sys.argv[0]} version:{__version__}",
        "date": str(datetime.datetime.now()),
        "reference": reffna,
        "reference_checksum": md5(reffna),
        "reference_features": [featureid for featureid in seq_d.keys()],
        "data": [{featureid: count_d[featureid] for featureid in count_d}],
    }
    data_d["data"][0]["id"] = sampleid
    print("data_d created.")

    with gzip.open("bam2counts.counts.json.gz", "wt") as fo:
        json.dump(data_d, fo, indent=4)
    print("json file created.")
