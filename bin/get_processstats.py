#!/usr/bin/env python3

import os
import sys
import datetime
import json
import gzip
import locale

__author__ = "Matthieu Pichaud"
__version__ = "0.3"

def parse_host(fn):

    # reads_splitmapped.py    input           2.qc/ERR275252_qced_pair_1.fq.gz        20729
    # reads_splitmapped.py    input           2.qc/ERR275252_qced_pair_2.fq.gz        20729
    # reads_splitmapped.py    output  inref   3.host/ERR275252_host_pair_inref_1.fq.gz        2980
    # reads_splitmapped.py    output  notinref        3.host/ERR275252_host_pair_notinref_1.fq.gz     17749
    # reads_splitmapped.py    output  inref   3.host/ERR275252_host_pair_inref_2.fq.gz        2980
    # reads_splitmapped.py    output  notinref        3.host/ERR275252_host_pair_notinref_2.fq.gz     17749
    # reads_splitmapped.py    input           2.qc/ERR275252_qced_single_1.fq.gz      1607
    # reads_splitmapped.py    output  inref   3.host/ERR275252_host_single1_inref_1.fq.gz     0
    # reads_splitmapped.py    output  notinref        3.host/ERR275252_host_single1_notinref_1.fq.gz  1607
    # reads_splitmapped.py    input           2.qc/ERR275252_qced_single_2.fq.gz      1569
    # reads_splitmapped.py    output  inref   3.host/ERR275252_host_single2_inref_1.fq.gz     1
    # reads_splitmapped.py    output  notinref        3.host/ERR275252_host_single2_notinref_1.fq.gz  1568

    stats_d = {
        "host_input_pair": [],
        "host_output_inref_pair": [],
        "host_output_notinref_pair": [],
        "host_input_single": [],
        "host_output_inref_single": [],
        "host_output_notinref_single": [],
    }

    with open(fn) as fi:
        for line in fi:
            items = line.rstrip("\r\n").split()
            inout = items[1]
            if "single" in items[3]:
                readtype = "single"
            elif (
                "pair_1" in items[3]
                or "pair_inref_1" in items[3]
                or "pair_notinref_1" in items[3]
            ):
                readtype = "pair"
            else:
                continue
            if inout == "input":
                key = f"host_{inout}_{readtype}"
            else:
                key = f"host_{inout}_{items[2]}_{readtype}"
            stats_d[key].append(int(items[-1]))

    for key in stats_d.keys():
        stats_d[key] = sum(stats_d[key])

    return stats_d

#
#     trigger_dest_d = {
#         "reads_splitmapped.py\tinput\t\t": ["host_input", -1],
#         "reads_splitmapped.py\toutput\tinref\t": ["host_inref", -1],
#         "reads_splitmapped.py\toutput\tnotintref\t": ["host_notinref", -1],
#     }
#
#     stats_d = {}
#     with open(fn) as fi:
#         for line in fi:
#             for trigger in trigger_dest_d.keys():
#                 if line.strip().startswith(trigger):
#                     dest = trigger_dest_d[trigger]
#                     items = line.rstrip("\r\n").split()
#                     stats_d[dest[0]] = stats_d.get(dest[0], []) + [
#                         int(items[dest[1]].replace(",", ""))
#                     ]
#                     break
#     return stats_d

def parse_trim_galore(fn):
    trigger_dest_d = {
        "Total number of sequences analysed:": ["trim_galore_input", -1],
        "Reads with adapters:": ["trim_galore_adapter", -2],
        "Number of sequence pairs removed because at least one read was shorter than the length cutoff": [
            "trim_galore_tooshort",
            -2,
        ],
    }

    stats_d = {}
    with open(fn) as fi:
        for line in fi:
            for trigger in trigger_dest_d.keys():
                if line.strip().startswith(trigger):
                    dest = trigger_dest_d[trigger]
                    items = line.rstrip("\r\n").split()
                    stats_d[dest[0]] = stats_d.get(dest[0], []) + [
                        int(items[dest[1]].replace(",", ""))
                    ]
                    break
    return stats_d

def parse_trimmomatic(fn):
    # Input Read Pairs: 24961 Both Surviving: 20729 (83.05%) Forward Only Surviving: 1607 (6.44%) Reverse Only Surviving: 1569 (6.29%) Dropped: 1056 (4.23%)
    stats_d = {}
    with open(fn) as fi:
        for line in fi:
            if line.startswith("Input Read Pairs:"):
                items = line.rstrip("\r\n").split()
                stats_d["trimmomatic_input"] = int(items[3])
                stats_d["trimmomatic_pairs"] = int(items[6])
                stats_d["trimmomatic_single"] = int(items[11]) + int(items[16])
                stats_d["trimmomatic_dropped"] = int(items[19])
                break
    return stats_d

def parse_kneaddata(fn):
    trigger_dest_d = {
        "Initial number of reads": "kneaddata_input",
        "Total reads after trimming": "kneaddata_posttrim",
        "Total reads after removing those found in reference database": "kneaddata_afterref",
        "Total reads after merging results from multiple databases": "kneaddata_aftermultiref",
    }
    stats_d = {}
    with open(fn) as fi:
        for line in fi:
            for trigger in trigger_dest_d.keys():
                if line.startswith(trigger):
                    dest = trigger_dest_d[trigger]
                    items = line.rstrip("\r\n").split()
                    stats_d[dest] = stats_d.get(dest, []) + [
                        [items[-3].split(os.sep)[-1], int(float(items[-1]))]
                    ]
                    break
    return stats_d

def parse_map(fn):
    stats_d = {}
    with open(fn, "rt") as fi:
        for iline, line in enumerate(fi):
            if iline == 0:
                stats_d["map_Npair"] = [int(line.rstrip())]
            elif iline == 1:
                stats_d["map_Nsingle"] = [int(line.rstrip())]
    return stats_d

if __name__ == "__main__":

    print(sys.argv)

    sampleid = sys.argv[1]
    trim_galorelog_fn = sys.argv[2]
    trimmomaticlog_fn = sys.argv[3]
    host_fn = sys.argv[4]
    maplog_fn = sys.argv[5]

    trim_galore_stats_d = parse_trim_galore(trim_galorelog_fn)
    trimmomatic_stats_d = parse_trimmomatic(trimmomaticlog_fn)
    host_stats_d = parse_host(host_fn)
    maplog_stats_d = parse_map(maplog_fn)

    data_d = {
        "generated_by": f"{sys.argv[0]} version:{__version__}",
        "date": str(datetime.datetime.now()),
        "data": [{}],
    }
    for stats_d in [
        trim_galore_stats_d,
        trimmomatic_stats_d,
        host_stats_d,
        maplog_stats_d,
    ]:
        for k, v in stats_d.items():
            data_d["data"][0][k] = v
    data_d["data"][0]["id"] = sampleid
    print("data_d created.")

    # Export
    with gzip.open("stats.json.gz", "wt") as fo:
        json.dump(data_d, fo, indent=4)
    print("json file created.")
