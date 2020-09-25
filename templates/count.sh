#! /bin/bash

##
##  Script to extract counts from bam alignment file
##
##  Author:         Matthieu Pichaud
##  Date:           24/07/2019
##
##  Requirements:
##
##  Arguments:
##  - sampleid:     sample identifier
##  - bamfile:      alignment of the sequencing data on the reference
##  - reffna:       reference fna
##  - bam2jsoncounts.py:  path to script used to count reads
##
##  Output:
##  - ${sampleid}x${refname}.count.json.gz
##  - ${sampleid}x${refname}.neighbor.tsv.gz
##

set +e

threads=${task.cpus} 

#  - audit trail
echo "=========="
echo $userName
date
echo "=========="
bam2jsoncounts.py --version
echo "=========="

#   - get counts from bam files
bam2jsoncounts.py -t \${threads}  "${sampleid}" "${reffna}" ${bamfile}
mv bam2counts.counts.json.gz  "${sampleid}x${refname}.count.json.gz"
mv bam2counts.neighbor.tsv.gz "${sampleid}x${refname}.neighbor.tsv.gz"
