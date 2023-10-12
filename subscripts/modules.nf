process fastp {
    debug true
    conda 'bioconda::fastp=0.23.4'
    publishDir "${params.work_dir}/fastp/reports", pattern: "{*html,*json}", mode: "copy"
    publishDir "${params.work_dir}/fastp", pattern: "*fq", mode: "copy"
    input:
    tuple val(f_name), path(fa)
    output:
    tuple val(f_name), path(fa)
	tuple val(f_name), path("${f_name}_R{1,2}_trimmed.fq"), emit: trimmed
	path "*{json,html}", emit: reports
    script:
    """
    fastp \
    -i ${fa[0]} \
    -o ${f_name}_R1_trimmed.fq \
    -I ${fa[1]} \
    -O ${f_name}_R2_trimmed.fq \
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
    output:
    path "*"
    script:
    println "${ref.getName()} will be indexed as $ref_number"
    if ( params.bbmap_memory ) {
        """
        bbmap.sh \
        ref=${ref} \
        build=${ref_number} \
        path="${params.work_dir}/bbmap_index" \
        -eoom -Xmx${params.bbmap_memory}
        """
    } else {
        """
        bbmap.sh \
        ref=${ref} \
        build=${ref_number} \
        path="${params.work_dir}/bbmap_index" \
        """
    }
}


process bbmap {
    debug true
    conda "bioconda::bbmap=39.01 bioconda::samtools=1.17"
    cpus $params.bbmap_cpus
    publishDir "${params.work_dir}/${bbmap_dir}/", pattern: "*bam", mode: "copy"
    publishDir "${params.work_dir}/${bbmap_dir}/unmapped", pattern: "*unmapped.fa", mode: "copy"
    publishDir "${params.work_dir}/${bbmap_dir}/statistics", pattern: "*txt", mode: "copy"
    input:
    tuple val(name), path(reads)
    val ref
    val ref_number
    val bbmap_dir
    output:
    tuple val(name), path("*unconc*{1,2}*"), emit: bbmap_tuple
    path "*txt", emit: stats
    path "*"
    script:
    if ( params.bbmap_memory ) {
        """
        bbmap.sh \
        ref=${ref} \
        build=${ref_number} \
        path="${params.work_dir}/bbmap_index" \
        in=${reads[0]} \
        in2=${reads[1]} \
        out=${name}.bam \
        outu=${name}_R#_unmapped.fa \
        statsfile="${name}_stats.txt" \
        -eoom -Xmx${params.bbmap_memory}
        """
    } else {
        """
        bbmap.sh \
        ref=${ref} \
        build=${ref_number} \
        path="${params.work_dir}/bbmap_index" \
        in=${reads[0]} \
        in2=${reads[1]} \
        out=${name}.bam \
        outu=${name}_R#_unmapped.fa \
        statsfile="${name}_stats.txt" \
        """
    }
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
    conda "bioconda::megahit=1.2.9"
    publishDir "${params.work_dir}/megahit/${megahit_subdir}", mode: "copy"
    input:
    tuple val(name), path(reads)
    val megahit_subdir
    output:
    path "*"
    script:
    """
    megahit \
    -1 ${reads[0]} \
    -2 ${reads[1]} \
    -o . \
    -t $params.megahit_threads
    -m $params.megahit_memory
    """
}
