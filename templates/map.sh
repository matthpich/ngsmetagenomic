#! /bin/bash

##
##  Script to map paired-end + single fq.gz file on reference using bwa and filtering output
##
##  Author:         Matthieu Pichaud
##  Date:           01/07/2020
##
##  Requirements:
##  - bwa
##  - samtools
##  - pysamstats
##
##  Arguments:
##  - sampleid:     sample identifier
##  - pair1/pair2/single:  read1/read2/single fq.gz
##  - refbwaindex:  reference bwa index
##  - reffna:       reference fna
##  - refname:      reference name
##  - sam_filter:   path to scrupt used to filter alignments based on various criteria
##
##  Output:
##  - ${sampleid}x${refname}.bam
##  - ${sampleid}x${refname}.samstats.gz
##  - map.json
##  - .command.sh
##  - .command.out
##

set -e

threads=${task.cpus}

pip install --no-cache-dir pysamstats

#  - audit trail
echo "=========="
echo $userName
date
echo "=========="
bwa 2>&1 >/dev/null | head -n 4 |tail -n +2
echo \$md5val
echo "=========="
samtools --version
echo "=========="
sam_filter.py --version
echo "=========="

mkdir ref
mv *.{amb,ann,bwt,pac,sa} ref/

#   - map
echo "${sampleid}x${refname}.paired.bam"
Ninputpairs1=`read_count.py ${pair1}` 
#Ninputpairs2=`read_count.py ${pair2}` # if Ninputpairs1 != Ninputpairs2, bwa will complain anyways
bwa mem -t\${threads} "ref/${refbwaallindex_base}" ${pair1} ${pair2} | samtools view -Sh -F4 - |sam_filter.py -l 50 -c 0 -p .95 - |grep .|grep -v "^#" > tmp
cat tmp | samtools view -H  - > tmp.header.sam
cat tmp | samtools view -bS - > tmp.paired.prelim.bam
Nmappedpairs=`samtools view tmp.paired.prelim.bam | awk '\$2!="*"' | cut -f 1 | sort | uniq | wc -l`

echo "${sampleid}x${refname}.single.bam"
Ninputsingle=`read_count.py ${single}`
bwa mem -t\${threads} "ref/${refbwaallindex_base}" ${single}        | samtools view -Sh -F4 - |sam_filter.py -l 50 -c 0 -p .95 - |grep .|grep -v "^#" > tmp
cat tmp | samtools view -H  - > tmp.header.sam
cat tmp | samtools view -bS - > tmp.single.prelim.bam
Nmappedsingle=`samtools view tmp.single.prelim.bam | awk '\$2!="*"' | cut -f 1 | sort | uniq | wc -l`

samtools merge tmp.all.prelim.bam tmp.paired.prelim.bam tmp.single.prelim.bam
samtools reheader tmp.header.sam tmp.all.prelim.bam |samtools sort - > "${sampleid}x${refname}.bam"
samtools index "${sampleid}x${refname}.bam"

#   - learn from bam
pysamstats -f "${reffna}" -t variation "${sampleid}x${refname}.bam" |pigz -p \${threads} -c > "${sampleid}x${refname}.samstats.gz"

#   - report
echo "sampleid\tactivity\tref\tin/out\tfile(s)\tsingle/pair\tNreads_mapped" >> map.log
echo "${sampleid}\tmap.sh\t${refname}\tin\t${pair1},${pair2}\tpair\t\${Ninputpairs1}" >> map.log
echo "${sampleid}\tmap.sh\t${refname}\tout\t${sampleid}x${refname}.bam\tpair\t\${Nmappedpairs}" >> map.log
echo "${sampleid}\tmap.sh\t${refname}\tin\t${single}\tsingle\t\${Ninputsingle}" >> map.log
echo "${sampleid}\tmap.sh\t${refname}\tout\t${sampleid}x${refname}.bam\tsingle\t\${Nmappedsingle}" >> map.log

#   - clean up
rm -f tmp*
