#! /usr/bin/env python3
# -*- coding: ascii -*-

"""
Map on reference and filter reads that map.
If a read of a pair map, both read of the pair are discarded.

Usage:
    reads_filtermapped.py [-h] [-v] [--sampleid SAMPLEID] [--withheader] [--threads THREADS] [--outprefix OUTPREFIX] BWAINDEX FASTQ ...

Arguments:
  BWAINDEX                             reference bwa index
  FASTQ                                1 or 2 fastq files

Options:
  -h, --help                           show this help
  -v, --version                        show version
  --withheader                         write header
  -s SAMPLEID, --sampleid SAMPLEID     sample identifier [default: ""]
  -t THREADS, --threads THREADS        number of threads used [default: 1]
  -o OUTPREFIX, --outprefix OUTPREFIX  prefix to use for outputs [default: reads_splitmapped]
"""

__author__ = "Matthieu Pichaud"
__version__ = "0.3"

from docopt import docopt
import os
import sys

def flexopen(fn, mode="rt"):
    if fn.endswith(".gzip") or fn.endswith(".gz"):
        return gzip.open(fn, mode)
    else:
        return open(fn, mode)

def read_count(fn):
    with flexopen(fn) as fi:
        for i, line in enumerate(fi): pass
    Nreads = float(i+1)/4
    if str(Nreads).endswith(".0"): return int(Nreads)
    else: return(Nreads)


def filterreads_basedonbwa(fastq_l, bwaindex, outprefix, threads, sampleid="", printheader=0):

     ## Run bwa
    reads = " ".join(fastq_l)
    os.system(f"bwa mem -t{threads} {bwaindex} {reads} > tmp.sam 2> /dev/null")

    ## List mapped reads from bwa output
    tofilter = set()
    with open("tmp.sam", "rt") as fi:
        for line in fi:
            if line[0] == "@":
                continue
            items = line.rstrip("\r\n").split("\t")
            if items[2] != "*":
                tofilter.add(items[0])

    ## Filter reads
    if len(fastq_l) == 1:
        fileset_l = [
            [
                fastq_l[0],
                f"{outprefix}_inref_1.fq",
                f"{outprefix}_notinref_1.fq",
            ]
        ]
    elif len(fastq_l) == 2:
        fileset_l = [
            [
                fastq_l[0],
                f"{outprefix}_inref_1.fq",
                f"{outprefix}_notinref_1.fq",
            ],
            [
                fastq_l[1],
                f"{outprefix}_inref_2.fq",
                f"{outprefix}_notinref_2.fq",
            ],
        ]

    for fileset in fileset_l:

        ## Count number of reads in input
        Nlreads_input = read_count(fileset[0])

        ## Split input read file in mapped and unmapped reads
        with open(fileset[0], "rt") as fi,\
             open(fileset[1], "wt") as fo_in,\
             open(fileset[2], "wt") as fo_notin:

            line_l = []
            Nreads_inref = 0
            Nreads_notinref = 0

            for line in fi:
                line = line.rstrip("\r\n")
                line_l.append(line)
                if len(line_l) == 4:
                    if line_l[0][1:].split()[0] in tofilter:
                        print("\n".join(line_l), file=fo_in)
                        Nreads_inref += 1 
                    else:
                        print("\n".join(line_l), file=fo_notin)
                        Nreads_notinref += 1 
                    line_l = []
 
    if printheader:
        print("sampleid\tactivity\tref\tin/out\tfile(s)\tsingle/pair\tNreads_filtered")
    
    print("\t".join([
        sampleid,
        "reads_splitmapped.py",
        bwaindex,
        "in",
        ",".join([fs[0] for fs in fileset_l]),
        "single "if len(fileset_l) == 1 else "pair",
        str(Nreads_notinref)]))
    
    print("\t".join([
        sampleid,
        "reads_splitmapped.py",
        bwaindex,
        "out",
        ",".join([fs[2] for fs in fileset_l]),
        "single "if len(fileset_l) == 1 else "pair",
        str(Nreads_notinref)]))

    return [fn for fileset in fileset_l for fn in fileset[1:]]

##
##      Main
##

if __name__ == "__main__":

    arguments = docopt(
        __doc__, version="{} version:{}".format(os.path.basename(__file__), __version__)
    )

    threads = int(arguments["--threads"])
    withheader = arguments["--withheader"]
    sampleid = arguments["--sampleid"]
    outprefix = arguments["--outprefix"]
    bwaindex = arguments["BWAINDEX"]
    fastq_l = arguments["FASTQ"]

    if len(fastq_l) not in [1, 2]:
        raise ValueError("Number of reads file provided different from 1 or 2")

    filterreads_basedonbwa(fastq_l, bwaindex, outprefix, threads, sampleid=sampleid, printheader=withheader)