#! /bin/bash

##
##  Collect information from fastp and mapping 
##
##  Author:         Matthieu Pichaud
##  Date:           01/07/2020
##
##  Arguments:
##  - sampleid:     sample identifier
##  - refname
##  - fastp.json
##  - map.log
##
##  Output:
##  - ${sampleid}x${refname}.stats.json.gz
##  - .command.sh
##  - .command.out
##

set -e

#  - audit trail
echo "=========="
echo $USER
date
echo "=========="
echo "=========="

#get_processstats.py ${sampleid} trim_galore.log trimmomatic.log host.log map.log
#mv stats.json.gz ${sampleid}x${refname}.seqstats.json.gz
touch ${sampleid}x${refname}.seqstats.json
pigz ${sampleid}x${refname}.seqstats.json
