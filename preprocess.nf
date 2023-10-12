include {fastp} from "${projectDir}/subscripts/modules.nf"
include {multiqc} from "${projectDir}/subscripts/modules.nf"
include {fastqc as fqc_before; fastqc as fqc_after} from "${projectDir}/subscripts/modules.nf"


workflow preprocess {
    take:
    reads_ch 
    
    main:
    fqc_before( reads_ch, 'qc_before' )
    fastp( reads_ch )
    fqc_after( fastp.out.trimmed, 'qc_after' )

    //merge reports from multiple processes
    reports = fqc_after.out.reports.concat(fqc_before.out.reports, fastp.out.reports).collect()
    multiqc( reports, 'multiqc' )

    emit:
    fastp.out.trimmed
}


workflow {
    Channel
        .fromFilePairs( params.reads, checkIfExists: true )
        .set { reads_ch }
    preprocess( reads_ch )
}