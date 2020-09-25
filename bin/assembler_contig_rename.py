#! /usr/bin/env python3

import string
import sys

# >NODE_1_length_580159_cov_136.902057

assembler_contig_f = sys.argv[1]
prefix = sys.argv[2]
assembler_contig_renamed_f = sys.argv[3]
assembler_contig_info_f = sys.argv[4]

with open(assembler_contig_f) as fi:
    line = fi.readline()
    if "multi=" in line and "len=" in line and "flag=" in line:
        assembler_source = "megahit"
    else:
        assembler_source = "spades"

with open(assembler_contig_f, "rt") as fi, open(
    assembler_contig_renamed_f, "wt"
) as fo1, open(assembler_contig_info_f, "wt") as fo2:

    if assembler_source == "spades":
        print("contig\tlength\tcoverage", file=fo2)
        for line in fi:
            line = line.rstrip("\r\n")
            if line[0] == ">":
                items = line[1:].split("_")
                contig_i = items[1]
                contig_len = items[3]
                contig_cov = items[5]
                print(">{}_{}".format(prefix, contig_i), file=fo1)
                print(
                    "{}_{}\t{}\t{}".format(prefix, contig_i, contig_len, contig_cov),
                    file=fo2,
                )
            else:
                print(line, file=fo1)

    elif assembler_source == "megahit":
        print("contig\tlength", file=fo2)
        for line in fi:
            line = line.rstrip("\r\n")
            if line[0] == ">":
                items = line[1:].rstrip("\r\n").split()
                contig_i = items[0].split("_")[1]
                contig_len = items[3][4:]
                contig_cov = "NA"
                print(">{}_{}".format(prefix, contig_i), file=fo1)
                print("{}_{}\t{}".format(prefix, contig_i, contig_len), file=fo2)
            else:
                print(line, file=fo1)
