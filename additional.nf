include {metawrap_binrefinement} from "${projectDir}/subscripts/process.nf"
include {gtdbtk_classify} from "${projectDir}/subscripts/process.nf"
include {checkm} from "${projectDir}/subscripts/process.nf"


workflow {
    if ( params.concoct ) {
        Channel
                .fromPath( "${workflow.launchDir}/${params.output}/metawrap_bins/concoct_bins/bin*" )
                .ifEmpty( "null1")
                .set {concoct_bins}
    } else {
        Channel 
                .fromPath( 'null1' )
                .set {concoct_bins}
    }
    if ( params.maxbin2 ) {
        Channel
                .fromPath( "${workflow.launchDir}/${params.output}/metawrap_bins/maxbin2_bins/bin*" )
                .ifEmpty( "null1" )
                .set {maxbin2_bins}
    } else {
        Channel 
                .fromPath( 'null2' )
                .set {maxbin2_bins}
    }
    if ( params.metabat2 ) {
        Channel
                .fromPath( "${workflow.launchDir}/${params.output}/metawrap_bins/metabat2_bins/bin*" )
                .ifEmpty( "null1" )
                .set {metabat2_bins}
    } else {
        Channel 
                .fromPath( 'null3' )
                .set {metabat2_bins}
    }
    
    metawrap_binrefinement( metabat2_bins.collect(), concoct_bins.collect(), maxbin2_bins.collect() )
    gtdbtk_classify( metawrap_binrefinement.out.refined_bins )
}
