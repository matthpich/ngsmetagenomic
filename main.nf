#!/usr/bin/env nextflow


// =============================================================================================
//     Parameters
// =============================================================================================

params.query            = ""
params.seqname_prefix   = ""
params.seqname_suffix   = ""
params.seqname_middle   = ""
params.seqNmax          = 9999

params.fastp_arguments  = "--qualified_quality_phred 15 --unqualified_percent_limit 40 --length_required 50 --low_complexity_filter 30 --cut_right --cut_right_window_size 4 --cut_right_mean_quality 15"
params.hostbwaindex     = ""
params.refbwaindex      = ""
params.reffna           = ""
params.refname          = "IGC"
params.jsoncount_derivative_ref = ""
params.richness_target  = 100000
params.richness_iter    = 30
params.MSP_descr        = ""
params.outfolder        = ""


// =============================================================================================
//     Facilitate access to parameters
// =============================================================================================

query                   = params.query
seqname_prefix          = params.seqname_prefix
seqname_suffix          = params.seqname_suffix
seqname_middle          = params.seqname_middle
seqNmax                 = params.seqNmax

fastp_arguments         = params.fastp_arguments
hostbwaindex            = params.hostbwaindex
refbwaindex             = params.refbwaindex
reffna                  = file(params.reffna)
refname                 = params.refname
jsoncount_derivative_ref = file(params.jsoncount_derivative_ref)
richness_target         = params.richness_target
richness_iter           = params.richness_iter
MSP_descr               = file(params.MSP_descr)
outfolder               = params.outfolder

userName                = workflow.userName


// =============================================================================================
//     Channels
// =============================================================================================

// Check if query a folder containing sequencing data or a SRA id ?

File query_path = new File(query)
println(query)
println(query_path)
println(query_path.isDirectory())

switch (query) {

    case { query_path.isDirectory() or query.startsWith("s3://") or query.startsWith("S3://") }:
        seqdataChannel = Channel
            .fromFilePairs("${query_path}/${seqname_prefix}*${seqname_middle}{1,2}${seqname_suffix}", flat: false)
            .take(seqNmax)
        break

    case { query.startsWith("SRP") or query.startsWith("ERR") or query.startsWith("ERP")}:
        seqdataChannel = Channel
            .fromSRA(query)
            .take(seqNmax)
        break

    default:
        println("I don't know");
        break
}

// Create reference index channels

hostbwaindexChannel     = Channel
                            .fromPath("${hostbwaindex}/*")

refbwaindexChannel      = Channel
                            .fromPath("${refbwaindex}/*")

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
    label       "qc"
    //storeDir  "${outfolder}/1.qc/${sampleid}"
    publishDir  "${outfolder}/1.qc/${sampleid}", mode: 'copy'

    input:
        set val(sampleid), file(reads_l) from seqdataChannel
        file hostbwaallindex from hostbwaindexChannel.collect()

    output:
        set val(sampleid), file("${sampleid}_QCed_pair1.fq.gz"), file("${sampleid}_QCed_pair2.fq.gz"), file("${sampleid}_QCed_single.fq.gz") into qceddataChannel
        set val(sampleid), file("fastp.json"), file("host*.log") into seqstatsChannel1
        file(".command.sh")
        //file(".command.log")
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


process map {

    tag       { "${sampleid}" }
    label     "map"
    //storeDir "${outfolder}/2.mapx${refname}/${sampleid}"
    publishDir "${outfolder}/2.mapx${refname}/${sampleid}", mode: 'copy'

    input:
        set val(sampleid), file(pair1), file(pair2), file(single) from qceddataChannel
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
    label     "count"
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
    label     "byproducts"
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

/*
seqstatsChannel = seqstatsChannel1.join( seqstatsChannel2 )

process seqstats {

    tag       { "${sampleid}" }
    label     "seqstats"
    storeDir  "${outfolder}/5.seqstatsx${refname}/${sampleid}"

    input:
        set val(sampleid), file("*.json") from seqstatsChannel

    output:
        file("${sampleid}x${refname}.seqstats.json.gz")

    script:
        template 'collect_mapstats.sh'

}
*/
