#! /bin/bash

##
##  Script to compute derivatives from json count file
##
##  Author:         Matthieu Pichaud
##  Date:           24/07/2019
##
##  Requirements:
##
##  Arguments:
##  - sampleid:     sample identifier
##  - refname:      alignment of the sequencing data on the reference
##  - jsoncount_derivative_ref: file listing path to various derivatives sqlite3 tables, format: derivative_name db_filepathv\tdb_table\tderivative_columns
##  - MSP_descr:    MSP description file, format: msp_name\tmodule_name\talpha\tprevalence\tNxokyok Nxdetecty0
##  - richness_target: downsizing target used before genes are counted
##  - richness_iter: number of iterations of gene richness computation to be made
##  - jsoncount_derivatives.py script used to compute derivatives from json file
##  - jsoncount_generichness.py script used to compute gene richness from json file
##  - jsoncount_colsum.py script used to compute colsum (and a simple version of gene richness)
##  - jsonMSP_modulepresenceabsence.py script used to compute module presence/absence based on MSP descriptions
##
##  Output:
##  - [sampleid]x[refname].[derivative_name].json.gz
##  - [sampleid]x[refname].generichness.json.gz
##  - [sampleid]x[refname].colsum.json.gz
##  - [sampleid]x[refname].MSP.modulepresenceabsence.json.gz
##

ls -lh

#  - audit trail
echo "=========="
echo $USER
date
echo "=========="
echo ${jsoncount_derivative_ref}
echo ${richness_target} ${richness_iter}
echo "=========="
jsoncount_derivatives.py --version
jsoncount_generichness.py --version
jsoncount_colsum.py --version
jsonMSP_modulepresenceabsence.py --version
echo "=========="

#   - create derivatives
echo "====="
jsoncount_derivatives.py      ${sampleid}x${refname}.count.json.gz ${jsoncount_derivative_ref}           ${sampleid}x${refname}
echo "====="
jsoncount_generichness.py    ${sampleid}x${refname}.count.json.gz ${richness_target} ${richness_iter} ${sampleid}x${refname}
echo "====="
jsoncount_colsum.py          ${sampleid}x${refname}.count.json.gz -n ${richness_target}                 ${sampleid}x${refname}
echo "====="
jsonMSP_modulepresenceabsence.py ${sampleid}x${refname}.MSP.json.gz ${MSP_descr}                        ${sampleid}x${refname}
echo "====="

