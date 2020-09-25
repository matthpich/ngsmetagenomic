#! /bin/bash

##
##  Script for paired-end / single-end fastq.gz  QC
##
##  Author:         Matthieu Pichaud
##  Date:           30/06/2020
##
##  Requirements:
##  - fastp
##  - bwa
##
##  Arguments:
##  - sampleid:                         sample identifier
##  - reads_l:                          reads1/reads2 file
##  - fastp_arguments:                  arguments used by fastp
##  - hostbwaindex:                     path to host bwa index file
##
##  Output:
##  - ${sampleid}_QCed_pair1.fq.gz
##  - ${sampleid}_QCed_pair2.fq.gz
##  - ${sampleid}_QCed_single.fq.gz
##  - fastq.json
##  - host.log
##

set -e

threads=${task.cpus}

#  - audit trail

echo "=========="
echo $USER
date +%Y%m%d_%H%M%S
echo "=========="
echo "fastp version: "`fastp --version`
echo "fastp_arguments: "${fastp_arguments}
bwa 2>&1 >/dev/null | head -n 4 |tail -n +2
echo "=========="


#  - fastp QC

mkdir -p 1.fastp

if [ $Nreads == 2 ]
then

    fastp \
        --in1 ${reads1} --in2 ${reads2} \
        --out1 1.fastp/${sampleid}_qced_pair_1.fq --out2 1.fastp/${sampleid}_qced_pair_2.fq \
        --unpaired1 1.fastp/${sampleid}_qced_single_1.fq --unpaired2 1.fastp/${sampleid}_qced_single_2.fq \
        - w \${threads} \
        ${fastp_arguments}

    cat 1.fastp/${sampleid}_qced_single_{1,2}.fq > 1.fastp/${sampleid}_qced_single.fq

elif [ $Nreads == 1 ]
then

    fastp \
        --in1 ${reads1} \
        --out1 1.fastp/${sampleid}_qced_single.fq \
        - w \${threads} \
        ${fastp_arguments}

    touch 1.fastp/{sampleid}_qced_pair_{1,2}.fq
fi 

# - mapping on host genome

mkdir -p 2.host

reads_splitmapped.py -o 2.host/${sampleid}_host_pair   --sampleid ${sampleid} --withheader -t \${threads} ${hostbwaallindex_base} 1.fastp/${sampleid}_qced_pair_1.fq 1.fastp/${sampleid}_qced_pair_2.fq >> host.log
reads_splitmapped.py -o 2.host/${sampleid}_host_single --sampleid ${sampleid}              -t \${threads} ${hostbwaallindex_base} 1.fastp/${sampleid}_qced_single.fq >> host.log

# - prepare exit
mv 2.host/${sampleid}_host_pair_notinref_1.fq   ${sampleid}_QCed_pair1.fq 
mv 2.host/${sampleid}_host_pair_notinref_2.fq   ${sampleid}_QCed_pair2.fq 
mv 2.host/${sampleid}_host_single_notinref_1.fq ${sampleid}_QCed_single.fq
pigz -p \${threads} ${sampleid}_QCed_{pair1,pair2,single}.fq

# - clean up
rm -rf 1.fastp 2.host
rm -f ${hostbwaallindex_base}*

sleep 1
