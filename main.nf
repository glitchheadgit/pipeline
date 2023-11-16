include {preprocess} from "${projectDir}/01_preprocess.nf"
include {mapping} from "${projectDir}/03_mapping.nf"
include {classify_reads} from "${projectDir}/04_classify_reads.nf"
include {assembly_megahit} from "${projectDir}/05_assembly_megahit.nf"
include {binning_metawrap} from "${projectDir}/06_binning_metawrap.nf"
include {binrefinement_metawrap} from "${projectDir}/06_binrefine_metawrap.nf"
include {gtdbtk_classify} from "${projectDir}/07_classify_bins.nf"
include {gunzip} from "${projectDir}/subscripts/process.nf"


workflow {
    //preprocessing
    Channel
        .fromFilePairs( params.reads, checkIfExists: true )
        .set { reads_ch }
    preprocess( reads_ch )

    //mapping
    ref_list = params.ref.split( /,/ )
    ref_list = ref_list.collect{ new File(it).getCanonicalPath() }
    ref_number_list = params.ref_number.toString().split(/,/)
    mapping( preprocess.out.trimmed, ref_list, ref_number_list )

    //reads classification
    classify_reads( preprocess.out.trimmed )

    //megahit assembly
    preprocess.out.trimmed
    .multiMap {
        read1: it[1][0]
        read2: it[1][1]
    }
    .set { read1_read2_ch }
    assembly_megahit( read1_read2_ch.read1.collect(), read1_read2_ch.read2.collect() )
    gunzip(  read1_read2_ch.read1.collect(), read1_read2_ch.read2.collect()  )

    //metawrap binning
    binning_metawrap( assembly_megahit.out.assembly, gunzip.out.reads1, gunzip.out.reads2 )
    binning_metawrap.out.metabat2_bins.ifEmpty( 'null1' )
    binning_metawrap.out.concoct_bins.ifEmpty( 'null2' )
    binning_metawrap.out.maxbin2_bins.ifEmpty( 'null3' )
}
