#!/usr/bin/env python3
# -*- coding: ascii -*-

"""
Check that sam file is valid by performing the following tests:
- the percentage identity, alignment length and coverage match user-determine criteria
- the length of the query and quality strings are the same
- the length of the query and cigar are the same 

Usage:
  sam_filter.py [-h] [-v] [-l <mal>] [-c <mpc>] [-p <mpi>] (-|-i <input>)

Options:
  -h, --help
  -v, --version
  -i <sam>, --input=<sam>        Input.
  -l <mal>, --minalignlen=<mal>  Minimum alignment length [default: 45].
  -c <mpc>, --minperccov=<mpc>   Minimum percentage coverage [default: 0].
  -p <mpi>, --minpercid=<mpi>    Minimum percentage identity [default: .95].

Version history:
  - 1.1: modified Nids so that "N" nucleotides, which are identified as mismatches in MD entry, do not impact percid
  - 1.2: now relies on pysam
         still does not use hard-/soft-clip information

Sources of information
"""

__author__ = "Matthieu Pichaud"
__version__ = "1.2"

from docopt import docopt
import pysam
import os

def fai2len(fai_fn):
    seqlen_d = {}
    with open(fai_fn) as fi:
        for line in fi:
            items = line.rstrip("\r\n").split("\t")
            seqlen_d[items[0]] = int(items[1])
    return seqlen_d

def sam2filter(sam_fn, minalignlen=0, minperccov=0, minpercid=0, minmapq=5):

    ## Parse sam file
    with pysam.AlignmentFile(sam_fn, mode="rb") as fi:

        yield (fi.header)

        for i, aln in enumerate(fi):

            ## Control that there is indeed a hit to *something*
            if aln.reference_id == "*":
                continue

            ## Filter MAPQ < minmapq
            if aln.mapping_quality < minmapq:
                continue

            ## Control alignment length compared to read length, number of matches, percentage identity
            qlen = aln.query_length
            rlen = aln.reference_length
            cigar_d = {}
            for k, v in aln.cigartuples:
                cigar_d[k] = cigar_d.get(k, 0) + v
            Nmatches = int(cigar_d.get(0, 0))  # read from cigar M
            if Nmatches < minalignlen:
                print(f"# Nmatches {Nmatches}/{minalignlen} - {aln.tostring()}")
                continue
            percmatches = float(Nmatches) / rlen
            if percmatches < minperccov:
                print(f"# percmatches {percmatches}/{minperccov} - {aln.tostring()}")
                continue
            Nid = Nmatches - aln.get_tag("NM")
            percid = float(Nid) / rlen
            if percid < minpercid:
                print(f"# percid {percid}/{minpercid} - {aln.tostring()}")
                continue

            yield aln.tostring()

if __name__ == "__main__":

    ##
    ##   Parse arguments
    ##

    arguments = docopt(
        __doc__, version="{} version:{}".format(os.path.basename(__file__), __version__)
    )
    if arguments["-"]:
        fi = "-"
    else:
        fi = arguments["--input"]
    minalignlen = int(arguments["--minalignlen"])
    minperccov = float(arguments["--minperccov"])
    minpercid = float(arguments["--minpercid"])

    ##
    ##   Run
    ##

    for aln in sam2filter(
        fi, minalignlen=minalignlen, minperccov=minperccov, minpercid=minpercid
    ):
        print(aln)
