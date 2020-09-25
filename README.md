This nextflow workflow can be used to process shotgun metagenomic data and get gene counts and variations.


# Getting started with a test

- install nextflow, see https://www.nextflow.io

- install conda, see https://www.anaconda.com

- clone the pipeline in the ngsmetagenomic folder  
```
nextflow clone https://gitlab.com/matthpich/ngsmetagenomic.git ngsmetagenomic
```
- run the pipeline
```
nextflow run ngsmetagenomic/main.nf -params-file ngsmetagenomic/test/params_test.yaml -profile conda -resume
```
done!

## What happened?

The clone you downloaded contains:
- the nextflow workflow in the **main.nf** file
- a basic nextflow config file **nextflow.config**
- a **templates** folder that contains all the scripts called by nextflow
- a **bin** folder that contains various supporting scripts
- a **test** folder with:
    - 2 small paired -end sequencing files, in the folder **sequencing_data**
    - a dummy host reference and its bwa indexes, in the folder **dummy_host**
    - a dummy microbiome gene reference and its bwa indexes, in the folder **dummy_reference**
    - a **derivatives.test.tsv** that describes the count matrices we want to derive from te raw count matrix
    - a **params_test.yaml** that helps set various parameters for the program to run (they can be inserted in the command line but it may cumbersome when there are many folders)

The nextflow submitted the following tasks:
- qc with trim_galore for adaptor trimming (cutadapt) and read quality trimming (fastqc), trimmomatic for another round of read quality trimming and finally a filter on the reads that map the host reference (bwa and script)
- map the data on the microbial gene reference (bwa and custom script) and saved where mismatches occured (pysamstats)
- collect qc and map statistics (script)
- count the number of hits on the genes of the microbial gene reference (script)
- derive from the raw count matrix various derivatives such as total number of hits per sample, gene richness, counts at the MSP level, and presence/absence of accessory MSP modules

The output folder contains the outputs of the various tasks, with the commands run and outputs in **.command.sh**, **.command.log** and **.command.out**

Finally, the **work** folder that contains all intermediate files may/should be clean up 



# Prerequisites

- install nextflow, see https://www.nextflow.io/

- install conda (alternative option: docker)

- get host reference fna and bwa indexes
For instance, to prepare the filtering of human reads:
```
mkdir hg19
cd hg19
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
mkdir -p index/bwa 
bwa index -p index/bwa/host <(pigz -dc hg19.fa.gz)
cd -
```

- get microbiome gene reference fna and index and bwa indexes
For instance, to prepare the mapping on the Integrated Gene Catalogue (https://doi.org/10.1038/nbt.2942):
```
mkdir IGC
cd IGC
wget ftp://parrot.genomics.cn/gigadb/pub/10.5524/100001_101000/100064/1.GeneCatalogs/IGC.fa.gz
pigz -d IGC.fa.gz
mkdir -p index/bwa
bwa index -p index/bwa/IGC IGC.fa.gz
cd -
```

- get the MSP to gene table
For instance, to prepare the use of MSPs (https://doi.org/10.1093/bioinformatics/bty830):
```
# visit https://academic.oup.com/bioinformatics/article/35/9/1544/5106712#supplementary-data and download the Supplementary Data - zip file
cd IGC
mkdir -p annotation/MSP
unzip /Users/matthieupichaud/Downloads/bty830_supp.zip
mv "bty830_Supplementary Table 1 - IGC.1661_MSPs.tsv.zip" annotation/MSP/IGC.1661_MSPs.tsv.zip
unzip annotation/MSP/IGC.1661_MSPs.tsv.zip
cd -
```

- create table of derivatives
The table is a tab-delimited file with a header:
    - "derivative_name": name of the derivative
    - "db\_filepathv": path to the tab-delimited text file / sqlite file that contains the mapping between the gene\_name and the category
    - "db_table": name of the table (for sqlite db)
    - "derivative_name": the name(s) of the column(s) that contains the groups


- get the MSP statistics to compute presence / absence table

# To run your samples

Create a yaml parameter file for your samples:
- seqfolder: folder that contains the sequencing files formatted as "*_{1,2}.f*q.gz"
- trim\_galore arguments: arguments used by trim\_galore
- trimmomatic\_arguments: arguments used by trimmomatic
- hostbwaindex: folder that contains the bwa index for the host
- reffna: path of the reference fasta file
- refname: name of the reference
- jsoncount\_derivative\_ref: derivatives from the raw count matrix
- richness_target: downsize target
- richness_iter: number of repetition of the downsizing to perform
- MSP_descr: 
- outfolder:
- Nsamplesmax:


# Todo

[ ] join result tables  
[ ] aggregate/format process statistincs in multiQC  
[ ] explain outputs  
[ ] add requested ressource allocation in parameter file
[ ] track technical conditions


