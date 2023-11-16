include {gtdbtk_classify} from "${projectDir}/subscripts/process.nf"


workflow gtdbtk {
    take:
    refined_bins

    main:
    gtdbtk_classify( refined_bins )
}

workflow {
    Channel
        .fromPath( params.refined_bins )
        .collect()
        .set { refined_bins }
    
    gtdbtk_classify( refined_bins )
}
