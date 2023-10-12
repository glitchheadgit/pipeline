include {fastqc; multiqc} from  "${projectDir}/config/modules.nf"

workflow qc {
    take:
        data //tuple of val(filename), path(file)
        qc_dir
    main:
        fastqc( data, qc_dir )
        multiqc( fastqc.out.fqc_out.collect(), qc_dir )
}