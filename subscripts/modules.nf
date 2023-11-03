process fastp {
    debug true
    conda 'bioconda::fastp=0.23.4'
    publishDir "${params.work_dir}/fastp/reports", pattern: "{*html,*json}", mode: "copy"
    publishDir "${params.work_dir}/fastp", pattern: "*fq", mode: "copy"
    input:
    tuple val(f_name), path(fa)
    output:
    tuple val(f_name), path(fa)
	tuple val(f_name), path("${f_name}_R{1,2}_trimmed.fq.gz"), emit: trimmed
	path "*{json,html}", emit: reports
    script:
    """
    fastp \
    -i ${fa[0]} \
    -o ${f_name}_R1_trimmed.fq.gz \
    -I ${fa[1]} \
    -O ${f_name}_R2_trimmed.fq.gz \
    -z $params.compress_level \
    -V \
    -g \
    --poly_g_min_len $params.poly_g_min_len \
    -x \
    --poly_x_min_len $params.poly_x_min_len \
    -5 \
    -3 \
    -M $params.cut_mean_quality \
    -n $params.n_base_limit \
    -e $params.average_qual \
    -l $params.length_required \
    -c \
    -w $params.working_thread_n \
    -j ${f_name}_fastp.json \
    -h ${f_name}_fastp.html
    """
}


process fastqc {
    debug true
    conda 'bioconda::fastqc=0.12.1'
    publishDir "${params.work_dir}/${fqc_dir}", pattern: "*html", mode: "copy"
    publishDir "${params.work_dir}/${fqc_dir}/zip", pattern: "*zip", mode: "copy"
    input:
    tuple val(index), path(reads)
    val fqc_dir
    output:
    tuple val(index), path(reads)
    path "*", emit: reports
    script:
    """
    fastqc -t 6 *
    """
}


process multiqc {
    debug true
    conda 'bioconda::multiqc=1.16 conda-forge::imp=2.19.0'
    publishDir "${params.work_dir}/${mqc_dir}", mode: "copy"
    input:
    path(reports)
    val mqc_dir
    output:
    path "*"
    script:
    """
    multiqc --config ${projectDir}/config/multiqc.yaml *
    """
}


process bbmap_index {
    debug true
    conda "bioconda::bbmap=39.01"
    input:
    tuple path(ref), val(ref_number)
    script:
    println "${ref.getName()} will be indexed as $ref_number"
    if ( params.bbmap_memory ) {
        """
        bbmap.sh \
        ref=${ref} \
        build=${ref_number} \
        path="${workflow.launchDir}/${params.work_dir}" \
        -eoom -Xmx${params.bbmap_memory}
        """
    } else {
        """
        bbmap.sh \
        ref=${ref} \
        build=${ref_number} \
        path="${workflow.launchDir}/${params.work_dir}"
        """
    }
}


process bbmap {
    debug true
    maxForks 4
    conda "bioconda::bbmap=39.01 bioconda::samtools=1.18"
    cpus params.bbmap_cpus
    publishDir "${params.work_dir}/${bbmap_dir}/", pattern: "*bam", mode: "copy"
    publishDir "${params.work_dir}/${bbmap_dir}/unmapped", pattern: "*_unmapped.fq.gz", mode: "copy"
    publishDir "${params.work_dir}/${bbmap_dir}/mapped", pattern: "*_mapped.fq.gz", mode: "copy"
    publishDir "${params.work_dir}/${bbmap_dir}/statistics", pattern: "*txt", mode: "copy"
    input:
    tuple val(name), path(reads)
    val ref
    val ref_number
    val bbmap_dir
    output:
    tuple val(name), path("*R{1,2}_unmapped.fq.gz"), emit: bbmap_unmapped_tuple
    tuple val(name), path("*R{1,2}_mapped.fq.gz"), emit: bbmap_mapped_tuple
    path "*txt", emit: stats
    path "*"
    script:
    flags = ""
    if ( params.bbmap_memory ) {
        flags += "-eoom -Xmx${params.bbmap_memory} "
    }
    """
    bbmap.sh \
    ref=${ref} \
    build=${ref_number} \
    path="${workflow.launchDir}/${params.work_dir}" \
    in=${reads[0]} \
    in2=${reads[1]} \
    out=${name}.bam \
    outu=${name}_R#_unmapped.fq.gz \
    outm=${name}_R#_mapped.fq.gz \
    statsfile="${name}_stats.txt" \
    $flags
    """
}


process merge_bbmap_statistics {
    debug true
    conda "conda-forge::pandas=2.1.1"
    publishDir "${params.work_dir}/${bbmap_dir}/statistics", mode: "copy"
    input:
    path stats
    val bbmap_dir
    output:
    path "*"
    script:
    """
    python3 ${workflow.projectDir}/subscripts/merge_bbmap_stat.py *
    """
}


process spades {
    debug true
    maxForks 1
    conda "bioconda::spades=3.15.5"
    publishDir "${params.work_dir}/spades", mode: "copy"
    input:
    val index
    val reads
    output:
    path "*"
    script:
    groups = ""
    [reads, index].transpose().each {read, index -> groups += "--pe-1 $index ${read[0]} --pe-2 $index ${read[1]} "}
    if ( params.spades_mode ) {
        """
        spades.py \
        -o ${params.spades_mode} \
        --${params.spades_mode} \
        $groups
        """
    } else {
        """
        spades.py \
        -o default \
        $groups
        """
    }
}


process metaphlan { // metaphlan --install --bowtie2db <database folder> needed
    debug true
    conda "bioconda::metaphlan=4.06"
    memory params.metaphlan_memory
    publishDir "${params.work_dir}/metaphlan", mode: "copy"
    input:
    tuple val(name),path(reads)
    val db
    output:
    path "*", emit: metaphlan_output
    script:
    """
    metaphlan \
    ${reads[0]},${reads[1]} \
    --input_type fastq \
    --bowtie2out $db \
    -o ${name}.txt \
    -t $params.metaphlan_analysis_type
    --nproc $params.metaphlan_cpus
    """

}


process metaphlan_merge {
    debug true
    conda "bioconda::metaphlan=4.06"
    publishDir "${params.work_dir}/metaphlan", mode: "copy"
    input:
    path metaphlan_output
    output:
    path "merged_abundance_table.txt"
    script:
    """
    merge_metaphlan_tables.py *.txt > merged_abundance_table.txt
    """
}


process megahit {
    debug true
    maxForks 1
    conda "bioconda::megahit=1.2.9"
    publishDir "${params.work_dir}", mode: "copy"
    input:
    path reads1
    path reads2
    output:
    path "megahit/*"
    path "megahit/final.contigs.fa", emit: assembly
    script:
    reads1 = reads1.join(',')
    reads2 = reads2.join(',')
    """
    megahit \
    -1 ${reads1} \
    -2 ${reads2} \
    -o megahit \
    """
}


process metawrap_binning {
    debug true
    conda "${projectDir}/mw-env" //"ursky::metawrap-mg=1.3.2"
    maxForks 1
    publishDir "${params.work_dir}", mode: "copy"
    beforeScript """ export CHECKM_DATA_PATH="${projectDir}/checkm_data" """
    input:
    path assembly
    path megahit_reads1, stageAs: "*_1.fastq"
    path megahit_reads2, stageAs: "*_2.fastq"
    output:
    path "metawrap_bins/*"
    path "metawrap_bins/metabat2_bins", emit: metabat2_bins, optional: true
    path "metawrap_bins/concoct_bins", emit: concoct_bins, optional: true
    path "metawrap_bins/maxbin2_bins", emit: maxbin2_bins, optional: true
    script:
    """
    metawrap binning \
    --metabat2 --maxbin2 --concoct \
    --run-checkm \
    -a $assembly \
    -o metawrap_bins \
    $megahit_reads1 $megahit_reads2
    """
}
// --maxbin2 --concoct 

process metawrap_binrefinement {
    debug true
    conda "${projectDir}/mw-env"
    maxForks 1
    publishDir "${params.work_dir}", mode: "copy"
    beforeScript """ export CHECKM_DATA_PATH="${projectDir}/checkm_data" """
    input:
    path metabat2_bins
    path concoct_bins
    path maxbin2_bins
    output:
    path "metawrap_binref/*"
    path "metawrap_binref/bin_refinement/metawrap*bins/*fa", emit: refined_bins
    script:
    """
    metawrap bin_refinement \
    -o metawrap_binref \
    -A $metabat2_bins -B $concoct_bins -C $maxbin2_bins \
    -c $params.metawrap_completion -x $params.metawrap_contamination
    """
}


process gtdbtk_classify {
    conda "bioconda::gtdbtk=2.3.2"
    maxForks 1
    publishDir "${params.work_dir}", mode: "copy"
    beforeScript """ export GTDBTK_DATA_PATH="${projectDir}/gtdbtk_data" """
    input:
    path bin, stageAs: "bins/*"
    output:
    stdout
    path "gtdbtk_classify/*"
    script:
    """
    mkdir scratch
    gtdbtk classify_wf \
    --mash_db ${workflow.launchDir}/${params.work_dir}/mash_db \
    --genome_dir bins/ \
    --extension fa \
    --out_dir gtdbtk_classify \
    --scratch_dir scratch \
    --cpus 20
    """
}


process kaiju {
    debug true
    conda "bioconda::kaiju"
    maxForks 1
    publishDir "${params.work_dir}/kaiju", mode: "copy"
    input:
    path reads1
    path reads2
    output:
    path "*"
    script:
    reads1 = reads1.join(',')
    reads2 = reads2.join(',')
    """
    kaiju-multi \
    -z 5 \
    -t nodes.dmp \
    -f ${projectDir}/kaiju_data/kaiju_db*fmi \
    -i $reads1 \
    -j $reads2 \
    """
}
