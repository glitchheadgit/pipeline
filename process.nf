include {bbmap as bbmap_human; bbmap as bbmap_host} from "${projectDir}/subscripts/modules.nf"
include {merge_bbmap_statistics as merge_human; merge_bbmap_statistics as merge_host} from "${projectDir}/subscripts/modules.nf"


workflow process {
    take:
    trimmed_reads_ch 
    
    main:
    if ( params.ref_number.size() == 2 ) {
        bbmap_human( trimmed_reads_ch, params.ref[0], params.ref_number[0].toInteger(), params.human_dir )
        merge_human( bbmap_human.out.stats.collect(),params.human_dir )
        bbmap_host( trimmed_reads_ch, params.ref[1], params.ref_number[1].toInteger(), params.host_dir )
        merge_host( bbmap_host.out.stats.collect(), params.host_dir)
    } else if ( params.ref_number.size() == 1 ) {
        bbmap_host( trimmed_reads_ch, params.ref[0], params.ref_number[0].toInteger(), params.host_dir )
        merge_host( bbmap_host.out.stats.collect(), params.host_dir )
    } else {
        println "Wrong number of reference, specify it only for a host and a human(optional)."
    }
}


workflow {
    Channel
        .fromFilePairs( params.reads, checkIfExists: true )
        .set { reads_ch }
    process( reads_ch )
}