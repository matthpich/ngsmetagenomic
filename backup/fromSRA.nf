#!/usr/bin/env nextflow


// =============================================================================================
//     Parameters
// =============================================================================================

params.sra              = ""

params.trim_galore_arguments = "-phred33 --quality 0 --stringency 5 --length 50"
params.trimmomatic_arguments = "HEADCROP:15 SLIDINGWINDOW:4:15 MINLEN:50"
params.hostbwaindex     = "/"
params.refbwaindex      = "/db/dmp/microbiome/pub/IGC/index/bwa"
params.reffna           = "/db/dmp/microbiome/pub/IGC/IGC.fna"
params.refname          = "IGC"
params.jsoncount_derivative_ref = "/home/pichama3/repos/ngsmetagenomic/bin/IGC2derivatives.tsv"
params.richness_target  = 100000
params.richness_iter    = 30
params.MSP_descr        = "/home/pichama3/repos/ngsmetagenomic/bin/HMP2profilexMSP20180113.module_describe.gz"
params.outfolder        = ""

params.keepNlinesfastq  = 9999999999
params.Nsamplesmax      = 9999999999


// =============================================================================================
//     Facilitate access to parameters
// =============================================================================================

sra                     = params.sra

trim_galore_arguments   = params.trim_galore_arguments
trimmomatic_arguments   = params.trimmomatic_arguments
hostbwaindex            = params.hostbwaindex
refbwaindex             = params.refbwaindex
reffna                  = file(params.reffna)
refname                 = params.refname
jsoncount_derivative_ref = file(params.jsoncount_derivative_ref)
richness_target         = params.richness_target
richness_iter           = params.richness_iter
MSP_descr               = file(params.MSP_descr)
outfolder               = params.outfolder

keepNlinesfastq         = params.keepNlinesfastq
Nsamplesmax             = params.Nsamplesmax


// =============================================================================================
//     Channels
// =============================================================================================

seqdataChannel          = Channel.fromSRA("${sra}").take(Nsamplesmax)
hostbwaindexChannel     = Channel.fromPath("${params.hostbwaindex}/*.{amb,ann,bwt,pac,sa}")
refbwaindexChannel      = Channel.fromPath("${refbwaindex}/*.{amb,ann,bwt,pac,sa}")

Channel
        .fromPath(jsoncount_derivative_ref)
        .splitCsv(header:true, sep:"\t")
        .map{ row -> file(row.db_filepathv) }
        .set { derivativefileChannel }


// =============================================================================================
//     Processes
// =============================================================================================

process qc {

    tag         { "${sampleid}" }
    label       "multi_cpu"
    publishDir  "${outfolder}/1.qc/${sampleid}"

    input:
        set val(sampleid), file(reads_l) from seqdataChannel
        file hostbwaallindex from hostbwaindexChannel.collect()

    output:
        set val(sampleid), file("${sampleid}_QCed_pair1.fq.gz"), file("${sampleid}_QCed_pair2.fq.gz"), file("${sampleid}_QCed_single.fq.gz") into qceddataChannel
        set val(sampleid), file("trim_galore.log"), file("trimmomatic.log"), file("host.log") into seqstatsChannel1
        file(".command.sh")
        file(".command.log")
        file(".command.out")

    script:
        hostbwaallindex_base = hostbwaallindex[0].toString() - ~/.amb?/ - ~/.ann?/ - ~/.bwt?/ - ~/.pac?/ - ~/.sa?/
        Nreads = reads_l.collect().size() 
        if (Nreads == 2){
            reads1 = reads_l[0]
            reads2 = reads_l[1]
        } else if (Nreads == 1){
            reads1 = reads_l[0]
        }
        template 'qc.sh'

}

qceddataChannel.into { qceddataChannel1; qceddataChannel2 }
qceddataChannel1.subscribe {  println "$it"  }

process map {

    tag       { "${sampleid}" }
    label     "multi_cpu"
    publishDir "${outfolder}/2.mapx${refname}/${sampleid}"

    input:
        set val(sampleid), file(pair1), file(pair2), file(single) from qceddataChannel2
        file( refbwaallindex ) from refbwaindexChannel.collect()
        file( reffna )

    output:
        set val(sampleid), file("${sampleid}x${refname}.bam") into mapChannel
        set val(sampleid), file("map.log") into seqstatsChannel2
        file("${sampleid}x${refname}.samstats.gz")
        file(".command.sh")
        file(".command.out")

    script:
        refbwaallindex_base = refbwaallindex[0].toString() - ~/.amb?/ - ~/.ann?/ - ~/.bwt?/ - ~/.pac?/ - ~/.sa?/
        template 'map.sh'

}

process count {

    tag       { "${sampleid}" }
    label     "multi_cpu"
    storeDir "${outfolder}/3.countx${refname}/${sampleid}"

    input:
        set val(sampleid), file(bamfile) from mapChannel
        file( reffna )

    output:
        set val(sampleid), file("${sampleid}x${refname}.count.json.gz") into countjsonChannel
        file("${sampleid}x${refname}.neighbor.tsv.gz")

    script:
        template 'count.sh'

}

process byproducts {

    tag       { "${sampleid}" }
    //label     "some_mem"
    cpus 2
    memory "10 GB"
    storeDir  "${outfolder}/4.byproductsx${refname}/${sampleid}"

    input:
        set val(sampleid), file("${sampleid}x${refname}.count.json.gz") from countjsonChannel
	file(derivativefile) from derivativefileChannel.collect()
	file(jsoncount_derivative_ref)

    output:
        file("*.json.gz")
        file("*.json.fingerprint.gz")
        file(".command.sh")
        file(".command.out")

    script:
        template 'byproducts.sh'
}

seqstatsChannel = seqstatsChannel1.join( seqstatsChannel2 )

process seqstats {

    tag       { "${sampleid}" }
    label     "basic"
    storeDir  "${outfolder}/5.seqstatsx${refname}/${sampleid}"

    input:
        set val(sampleid), file("trim_galore.log"), file("trimmomatic.log"), file("host.log"), file("map.log") from seqstatsChannel

    output:
        file("${sampleid}x${refname}.seqstats.json.gz")

    script:
        template 'collect_mapstats.sh'

}
