include {kraken2} from "${projectDir}/subscripts/process.nf"
include {bracken} from "${projectDir}/subscripts/process.nf"
include {merge_kr2_br} from "${projectDir}/subscripts/process.nf"


workflow classify_reads {
    take:
    reads

    main:
    kraken2( reads )
    bracken( kraken2.out.kraken2_out )
    merge_kr2_br( kraken2.out.stats, bracken.out.stats )
}

workflow {
    Channel
        .fromFilePairs( params.reads, checkIfExists: true )
        .set { reads_ch }
    classify_reads( reads_ch )
}
