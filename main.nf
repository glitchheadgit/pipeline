include {preprocess} from "${projectDir}/submodules/01_preprocess.nf"
include {mapping} from "${projectDir}/submodules/03_mapping.nf"
include {classify_reads} from "${projectDir}/submodules/04_classify_reads.nf"
include {assembly_megahit} from "${projectDir}/submodules/05_assembly_megahit.nf"
include {binning_metawrap} from "${projectDir}/submodules/06_binning_metawrap.nf"
include {binrefinement_metawrap} from "${projectDir}/submodules/06_binrefine_metawrap.nf"
include {gtdbtk_classify} from "${projectDir}/submodules/07_classify_bins.nf"

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
    preprocess.out.trimmed
        .multiMap {
            read1: it[1][0]
            read2: it[1][1]
        }
        .set { read1_read2_ch }
    classify_reads( read1_read2_ch.read1.collect(), read1_read2_ch.read2.collect() )
    //megahit assembly
    assembly_megahit( classify_reads.out.reads1, classify_reads.out.reads2 )
    //metawrap binning
    binning_metawrap( assembly_megahit.out.assembly, assembly_megahit.out.reads1, assembly_megahit.out.reads2 )
    //metawrap binrefinement
    binrefinement_metawrap( binning_metawrap.out.metabat2_bins, binning_metawrap.out.concoct_bins, binning_metawrap.out.maxbin2_bins )
    //gtdbtk classify
    gtdbtk_classify( binrefinement_metawrap.out.refined_bins )
}