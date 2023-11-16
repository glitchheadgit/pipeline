process fastp {
    debug true
    conda 'bioconda::fastp=0.23.4'
    publishDir "${params.output}/fastp/reports", pattern: "{*html,*json}", mode: "copy"
    publishDir "${params.output}/fastp", pattern: "*gz", mode: "copy"
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
    conda 'bioconda::fastqc=0.12.1'
    publishDir "${params.output}/${fqc_dir}", pattern: "*html", mode: "copy"
    publishDir "${params.output}/${fqc_dir}/zip", pattern: "*zip", mode: "copy"
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
    conda 'bioconda::multiqc conda-forge::imp=2.19.0'
    publishDir "${params.output}/${mqc_dir}", mode: "copy"
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
        path="${workflow.launchDir}" \
        -eoom -Xmx${params.bbmap_memory}
        """
    } else {
        """
        bbmap.sh \
        ref=${ref} \
        build=${ref_number} \
        path="${workflow.launchDir}"
        """
    }
}


process bbmap {
    debug true
    maxForks 4
    conda "bioconda::bbmap=39.01 bioconda::samtools=1.18"
    cpus params.bbmap_cpus
    publishDir "${params.output}/${bbmap_dir}/", pattern: "*bam", mode: "copy"
    publishDir "${params.output}/${bbmap_dir}/unmapped", pattern: "*_unmapped.fq.gz", mode: "copy"
    publishDir "${params.output}/${bbmap_dir}/mapped", pattern: "*_mapped.fq.gz", mode: "copy"
    publishDir "${params.output}/${bbmap_dir}/statistics", pattern: "*txt", mode: "copy"
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
    path="${workflow.launchDir}" \
    in=${reads[0]} \
    in2=${reads[1]} \
    out=${name}.bam \
    outu=${name}_R#_unmapped.fq.gz \
    outm=${name}_R#_mapped.fq.gz \
    statsfile="${name}_stats.txt" \
    $flags
    """
}


process kraken2 {
	conda = 'bioconda::kraken2'
	maxForks 1
	publishDir "${params.output}/kraken2", pattern: "*kreport", mode: "copy"
	publishDir "${params.output}/kraken2/unclassified_reads", pattern: "*unclassified*", mode: "copy"
	input:
    tuple val(index), path(reads)
	output:
    tuple val(index), path("*.kreport"), emit: kraken2_out
	path "*"
	script:
	"""
	kraken2 --db ${projectDir}/kraken2_data --report ${index}.kreport --gzip-compressed --paired --unclassified-out ${index}_unclassified_#.fq ${reads[0]} ${reads[1]}
	"""
}


process bracken {
    conda = 'bioconda::bracken'
    publishDir "${params.output}/bracken", pattern: "*bracken", mode: 'copy'
    input:
    tuple val(index), path(kreport)
    val class_lvl
    output:
    tuple val(index), path("*.bracken"), emit: bracken_out
    path "*"
    script:
    """
    bracken -d $params.kraken_db -i $kreport -o ${index}.bracken -l $class_lvl
    """
}


process merge_bbmap_statistics {
    conda "conda-forge::pandas=2.1.1"
    publishDir "${params.output}/${bbmap_dir}/statistics", mode: "copy"
    input:
    path stats
    val bbmap_dir
    output:
    path "*"
    script:
    """
    python3 ${workflow.projectDir}/subscripts/merge_bbmap_stat.py -i *
    """
}


process spades {
    debug true
    maxForks 1
    conda "bioconda::spades=3.15.5"
    publishDir "${params.output}/spades", mode: "copy"
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
    publishDir "${params.output}/metaphlan", mode: "copy"
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
    publishDir "${params.output}/metaphlan", mode: "copy"
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
    publishDir "${params.output}", mode: "copy"
    input:
    path reads1
    path reads2
    output:
    path "megahit/final.contigs.fa", emit: assembly
    path "megahit"
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


process metawrapbin {
    conda "${projectDir}/mw-env" //"ursky::metawrap-mg=1.3.2"
    maxForks 1
    publishDir "${params.output}", mode: "copy"
    input:
    path assembly
    path reads1, stageAs: 'read*_1.fastq'
    path reads2, stageAs: 'read*_2.fastq'
    output:
    path "metawrap_bins/*"
    path "metawrap_bins/metabat2_bins/*", emit: metabat2_bins, optional: true
    path "metawrap_bins/concoct_bins/*", emit: concoct_bins, optional: true
    path "metawrap_bins/maxbin2_bins/*", emit: maxbin2_bins, optional: true
    script:
    binners = ""
    if (params.concoct) {
        binners += "--concoct "
    }
    if (params.maxbin2) {
        binners += "--maxbin2 "
    }
    if (params.metabat2) {
        binners += "--metabat2 "
    }
    """
    echo ${projectDir}/checkm_data | checkm data setRoot
    metawrap binning \
    $binners\
    --run-checkm \
    -a $assembly \
    -o metawrap_bins \
    $reads1 $reads2
    """
}
// --maxbin2 --concoct 

process metawrap_binrefinement {
    debug true
    conda "${projectDir}/mw-env"
    maxForks 1
    errorStrategy "ignore"
    publishDir "${params.output}", mode: "copy"
    input:
    path metabat2_bins, stageAs: "metabat2/*"
    path concoct_bins, stageAs: "concoct/*"
    path maxbin2_bins, stageAs: "maxbin2/*"
    output:
    path "metawrap_binref/*"
    path "metawrap_binref/metawrap*bins/*fa", emit: refined_bins
    script:
    not_nulls = []
    bin_dirs = ''
    if (!(metabat2_bins ==~ /.*null.*/)) {
        not_nulls.add('metabat2')
    }
    if (!(concoct_bins ==~ /.*null.*/)) {
        not_nulls.add('concoct')
    }
    if (!(maxbin2_bins ==~ /.*null.*/)) {
        not_nulls.add('maxbin2')
    }
    if (not_nulls.size == 3) {
        bin_dirs += "-A ${not_nulls[0]} -B ${not_nulls[1]} -C ${not_nulls[2]}"
    }
    if (not_nulls.size == 2) {
        bin_dirs += "-A ${not_nulls[0]} -B ${not_nulls[1]}"
    }
    if (not_nulls.size == 1) {
        bin_dirs += "-A ${not_nulls[0]}"
    }
    """
    metawrap bin_refinement \
    -o metawrap_binref \
    -t 8 \
    -m 100 \
    -c $params.metawrap_completion -x $params.metawrap_contamination \
    $bin_dirs
    """
}


process checkm {
    debug true
    conda "${projectDir}/mw-env"
    publishDir "metawrap_binref/bin_refinement/metawrap*bins/", mode: "copy"
    input:
    path refined_bins, stageAs: "metawrap_bins/*"
    output:
    path "*"
    script:
    """
    bash ${projectDir}/mw-env/bin/metawrap-scripts/run_checkm.sh metawrap_bins
    """
}


process gunzip {
    input:
    path "R1_*.fq.gz"
    path "R2_*.fq.gz"
    output:
    path "R1*.fq", emit: reads1
    path "R2*.fq", emit: reads2
    script:
    """
    gunzip -k --force R*.gz
    """
}
process gtdbtk_classify {
    debug true
    conda "bioconda::gtdbtk=2.3.2"
    maxForks 1
    publishDir "${params.output}", mode: "copy"
    input:
    path bin, stageAs: "bins/*"
    output:
    stdout
    path "gtdbtk_classify/*"
    script:
    """
    export GTDBTK_DATA_PATH=${projectDir}/gtdbtk_data
    mkdir scratch
    gtdbtk classify_wf \
    --mash_db ${workflow.launchDir}/${params.output}/mash_db \
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
    publishDir "${params.output}/kaiju", mode: "copy"
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
