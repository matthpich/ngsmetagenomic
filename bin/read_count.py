#!/usr/bin/env python3

import sys
import gzip


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

##
##      Main
##

if __name__ == "__main__":
    print(read_count(sys.argv[1]))