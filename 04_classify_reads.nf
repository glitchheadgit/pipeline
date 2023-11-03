include {kaiju} from "${projectDir}/subscripts/modules.nf"


workflow classify_reads {
    take:
    reads1
    reads2

    main:
    kaiju( reads1, reads2 )

    emit:
    reads1
    reads2
}

workflow {
    Channel
        .fromFilePairs( params.reads, checkIfExists: true )
        .multiMap {
            read1: it[1][0]
            read2: it[1][1]
        }
        .set { reads_ch }
    classify_reads( reads_ch.read1.collect(), reads_ch.read2.collect() )
}
