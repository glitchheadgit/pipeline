include {bbmap as bbmap_human; bbmap as bbmap_host} from "${projectDir}/subscripts/process.nf"
include {merge_bbmap_statistics as merge_statistics_human; merge_bbmap_statistics as merge_statistics_host} from "${projectDir}/subscripts/process.nf"


workflow mapping {
    take:
    trimmed_reads_ch 
    ref_list
    ref_number_list

    main:
    if ( ref_number_list.size() == 2 && ref_list.size() == 2 ) {
        bbmap_human( trimmed_reads_ch, ref_list[0], ref_number_list[0], params.human_dir )
        merge_statistics_human( bbmap_human.out.stats.collect(),params.human_dir )
        bbmap_host( trimmed_reads_ch, ref_list[1], ref_number_list[1], params.host_dir )
        merge_statistics_host( bbmap_host.out.stats.collect(), params.host_dir)
    } else if ( ref_number_list.size() == 1 && ref_list.size() == 1 ) {
        bbmap_host( trimmed_reads_ch, ref_list[0], ref_number_list[0], params.host_dir )
        merge_statistics_host( bbmap_host.out.stats.collect(), params.host_dir )
    } else {
        println "Wrong number of references or reference numbers, specify it only for a host and a human(optional)."
    }
}


workflow {
    ref_list = params.ref.split(/,/)
    ref_list = ref_list.collect{ new File(it).getCanonicalPath() }
    ref_number_list = params.ref_number.toString().split(/,/)
    Channel
        .fromFilePairs( params.reads, checkIfExists: true )
        .set { reads_ch }

    mapping( reads_ch, ref_list, ref_number_list )
}
