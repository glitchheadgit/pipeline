include {metawrap_binning} from "${projectDir}/subscripts/modules.nf"


workflow binning_metawrap {
    take:
    assembly
    reads1
    reads2

    main:
    metawrap_binning( assembly, reads1, reads2 )

    emit:
    metabat2_bins = metawrap_binning.out.metabat2_bins
    concoct_bins = metawrap_binning.out.concoct_bins
    maxbin2_bins = metawrap_binning.out.maxbin2_bins
}

workflow {
    Channel
        .fromFilePairs( params.reads, checkIfExists: true )
        .multiMap {
            read1: it[1][0]
            read2: it[1][1]
        }
        .set { reads_ch }
    assembly = new File( params.assembly )
    binning_metawrap( assembly, reads_ch.read1.collect(), reads_ch.read2.collect() )
}