#! /usr/bin/env python3

import sys
import os
import pysam
import gzip

prefix = sys.argv[1]
if len(sys.argv) == 2:
    f = sys.stdin
elif len(sys.argv) == 3:
    f = sys.argv[2]
cwd = os.getcwd()

with pysam.AlignmentFile(f, "r") as f, gzip.open(
    os.path.join(cwd, prefix + ".mapped.pair.1.fq.gz"), "wt"
) as fmp1, gzip.open(
    os.path.join(cwd, prefix + ".mapped.pair.2.fq.gz"), "wt"
) as fmp2, gzip.open(
    os.path.join(cwd, prefix + ".mapped.single.1.fq.gz"), "wt"
) as fms1, gzip.open(
    os.path.join(cwd, prefix + ".mapped.single.2.fq.gz"), "wt"
) as fms2, gzip.open(
    os.path.join(cwd, prefix + ".unmapped.pair.1.fq.gz"), "wt"
) as fup1, gzip.open(
    os.path.join(cwd, prefix + ".unmapped.pair.2.fq.gz"), "wt"
) as fup2, gzip.open(
    os.path.join(cwd, prefix + ".unmapped.single.1.fq.gz"), "wt"
) as fus1, gzip.open(
    os.path.join(cwd, prefix + ".unmapped.single.2.fq.gz"), "wt"
) as fus2:

    case2file_d = {
        (1, 1, 1): [fmp1, 0],
        (1, 1, 2): [fmp2, 0],
        (1, 0, 1): [fms1, 0],
        (1, 0, 2): [fms2, 0],
        (0, 1, 1): [fus1, 0],
        (0, 1, 2): [fus2, 0],
        (0, 0, 1): [fup1, 0],
        (0, 0, 2): [fup2, 0],
    }

    for aln in f:

        # read 1 or 2
        if aln.is_read1:
            iread = 1
        else:
            iread = 2
        # read mapped or not
        if aln.is_unmapped:
            ismap = 0
        else:
            ismap = 1
        # paired read is mapped or not
        if aln.mate_is_unmapped:
            ismatemap = 0
        else:
            ismatemap = 1

        print(
            "@%s\n%s\n+\n%s"
            % (
                aln.query_name,
                aln.query_sequence,
                "".join([chr(x + 33) for x in aln.query_qualities]),
            ),
            file=case2file_d[(ismap, ismatemap, iread)][0],
        )
        case2file_d[(ismap, ismatemap, iread)][1] += 1

## Report
print("prefix\tun/mapped\tsingle/pair\t1/2\tN")
print(f"{prefix}\tmapped\tpair\t1\t{case2file_d[1,1,1][1]}")
print(f"{prefix}\tmapped\tpair\t2\t{case2file_d[1,1,2][1]}")
print(f"{prefix}\tmapped\tsingle\t1\t{case2file_d[1,0,1][1]}")
print(f"{prefix}\tmapped\tsingle\t2\t{case2file_d[1,0,2][1]}")
print(f"{prefix}\tunmapped\tpair\t1\t{case2file_d[0,0,1][1]}")
print(f"{prefix}\tunmapped\tpair\t2\t{case2file_d[0,0,2][1]}")
print(f"{prefix}\tunmapped\tsingle\t1\t{case2file_d[0,1,1][1]}")
print(f"{prefix}\tunmapped\tsingle\t2\t{case2file_d[0,1,2][1]}")
