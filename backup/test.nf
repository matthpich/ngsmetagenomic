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

process qc {

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
        set val(sampleid), file("${sampleid}_QCed_pair1.fq.gz"), file("${sampleid}_QCed_pair2.fq.gz"), file("${sampleid}_QCed_single.fq.gz") into qceddataChannel
        set val(sampleid), file("trim_galore.log"), file("kneaddata.log") into seqstatsChannel1
        file(".command.sh")
        file(".command.log")
        file(".command.out")

    script:
        template 'qc.sh'

}

process map {

    tag       { "${sampleid}" }
    label     "multi_cpu"
    storeDir  "${outfolder}/2.mapx${refname}/${sampleid}"
    beforeScript "source /db/dmp/microbiome/pub/profiles/src/NGSpublicxmap.dev.setenv.sh"
    errorStrategy "retry"
    maxRetries 3

    input:
        set val(sampleid), file(pair1), file(pair2), file(single) from qceddataChannel

    output:
        set val(sampleid), file("${sampleid}x${refname}.bam") into mapChannel
        set val(sampleid), file("map.log") into seqstatsChannel2
        file("${sampleid}x${refname}.samstats.gz")
        file(".command.sh")
        file(".command.out")

    script:
        template 'map.sh'

}

process count {

    tag       { "${sampleid}" }
    label     "multi_cpu"
    storeDir "${outfolder}/3.countx${refname}/${sampleid}"
    beforeScript "source /db/dmp/microbiome/pub/profiles/src/NGSpublicxmap.dev.setenv.sh"
    errorStrategy "retry"
    maxRetries 3

    input:
        set val(sampleid), file(bamfile) from mapChannel

    output:
        set val(sampleid), file("${sampleid}x${refname}.count.json.gz") into countjsonChannel
        file("${sampleid}x${refname}.neighbor.tsv.gz")

    script:
        template 'count.sh'

}

process byproducts {

    tag       { "${sampleid}" }
    label     "basic"
    storeDir  "${outfolder}/4.byproductsx${refname}/${sampleid}"
    beforeScript "source /db/dmp/microbiome/pub/profiles/src/NGSpublicxmap.dev.setenv.sh"
    errorStrategy "retry"
    maxRetries 3

    input:
        set val(sampleid), file( "${sampleid}x${refname}.count.json.gz" ) from countjsonChannel

    output:
        file ( "*.json.gz" )
        file( ".command.sh" )
        file( ".command.out" )

    script:
        template 'byproducts.sh'
}

seqstatsChannel = seqstatsChannel1.join( seqstatsChannel2 )

process seqstats {

    tag       { "${sampleid}" }
    label     "basic"
    storeDir  "${outfolder}/5.seqstatsx${refname}/${sampleid}"
    beforeScript "source /db/dmp/microbiome/pub/profiles/src/NGSpublicxmap.dev.setenv.sh"
    errorStrategy "retry"
    maxRetries 3

    input:
        set val(sampleid), file("trim_galore.log"), file("kneaddata.log"), file("map.log") from seqstatsChannel

    output:
        file("${sampleid}x${refname}.seqstats.json.gz")
    
    script:
        template 'collect_mapstats.sh'                


}

