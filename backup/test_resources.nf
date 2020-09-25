#!/usr/bin/env nextflow

params.seqfolder        = "/home/pichama3/repos/ngsmetagenomic/data/"
params.hostbowtie2      = "/db/dmp/microbiome/ref/hg19/index/bowtie2/hg19"
params.refbwaindex      = "/db/dmp/microbiome/pub/IGC/index/bwa/IGC"
params.reffna           = "/db/dmp/microbiome/pub/IGC/IGC.fna"
params.refname          = "IGC"
params.jsoncount_derivative_ref = "/home/pichama3/repos/ngsmetagenomic/bin/IGC2derivatives.tsv"
params.richness_target  = 100000
params.richness_iter    = 30
params.MSP_descr        = "/home/pichama3/repos/ngsmetagenomic/bin/HMP2profilexMSP20180113.module_describe.gz"
params.outfolder        = "/home/pichama3/repos/ngsmetagenomic/output/"
params.threads          = 30

seqfolder               = file(params.seqfolder)
hostbowtie2             = file(params.hostbowtie2)
refbwaindex             = file(params.refbwaindex)
reffna                  = file(params.reffna)
refname                 = params.refname
jsoncount_derivative_ref = file(params.jsoncount_derivative_ref)
richness_target         = params.richness_target
richness_iter           = params.richness_iter
MSP_descr               = file(params.MSP_descr)
outfolder               = file(params.outfolder)
threads                 = params.threads

trim_galore_arguments   = "-phred33 --quality 0 --stringency 5 --length 50"
trimmomatic_arguments   = "HEADCROP:15 SLIDINGWINDOW:4:15 MINLEN:50"
hostDNA                 = file("/db/dmp/microbiome/ref/hg19/hg19")



seqdataChannel          = Channel
                            .fromFilePairs("${seqfolder}/*_{1,2}.f*q.gz", flat: true)

process qc_test {

    //conda     "python=3.6 bioconda::trim-galore=0.6.3 bioconda::kneaddata"
    tag       { "${sampleid}" }
    label     "multi_cpu"
    storeDir  "${outfolder}/1.qc/${sampleid}"
    beforeScript "source /db/dmp/microbiome/pub/profiles/src/NGSpublicxmap.dev.setenv.sh"
    errorStrategy "retry"
    maxRetries 3
 
    input:
        set val(sampleid), file(read1), file(read2) from seqdataChannel

    output:
        file("test.log")
        file(".command.sh")
        file(".command.log")
        file(".command.out")

    """
    samtools > test.log
    """
}

