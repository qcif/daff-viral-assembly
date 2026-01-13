//Original code
/*
process KRAKEN2 {
  tag "${sampleid}"
  label "setting_23"
  containerOptions "${bindOptions}"
  publishDir "${params.outdir}/${sampleid}/05_read_classification", mode: 'copy'

  input:
    tuple val(sampleid), path(fastq1), path(fastq2)
  output:
//    path("${sampleid}_kraken2.log")
//    path("${sampleid}_kraken2_report.txt")
//    path("${sampleid}_kraken2_output.txt")
//    path("${sampleid}_unclassified_1.fastq")
//    path("${sampleid}_unclassified_2.fastq")
    tuple val(sampleid), path("${sampleid}_kraken2_report.txt"), path("${sampleid}_kraken2_output.txt"), path(fastq1), path(fastq2), path("${sampleid}_unclassified_1.fastq"), path("${sampleid}_unclassified_2.fastq"), emit: kraken2_results
    tuple val(sampleid), path("${sampleid}_kraken2_report.txt"), emit: kraken2_results2
  script:
  """
  kraken2 --db ${params.kraken2_db} --use-names \
          --paired --threads ${task.cpus} \
          --gzip-compressed \
          --confidence 0.05 \
          --report ${sampleid}_kraken2_report.txt \
          --output ${sampleid}_kraken2_output.txt \
          --unclassified-out ${sampleid}_unclassified#.fastq \
          --report-minimizer-data \
          --minimum-hit-groups 3 \
          ${fastq1} ${fastq2} > ${sampleid}_kraken2.log
  
  """
}
*/
process KRAKEN2_KRAKEN2 {
    tag "$meta.id"
    label 'setting_23'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0f/0f827dcea51be6b5c32255167caa2dfb65607caecdc8b067abd6b71c267e2e82/data' :
        'community.wave.seqera.io/library/kraken2_coreutils_pigz:920ecc6b96e2ba71' }"

    input:
    tuple val(meta), path(reads)
    path db
    val save_classified_output_fastqs
    val save_unclassified_output_fastqs
    val save_reads_assignment

    output:
    tuple val(meta), path('*.classified{.,_}*')     , optional:true, emit: classified_reads_fastq
    tuple val(meta), path('*.unclassified{.,_}*')   , optional:true, emit: unclassified_reads_fastq
    tuple val(meta), path('*_kraken2_classified_reads.txt')   , optional:true, emit: classified_reads_assignment
    //tuple val(meta), path('*report.txt')                           , emit: report
    tuple val(meta.id), path('*report.txt')                           , emit: report
    path "versions.yml"                                            , emit: versions
    tuple val(meta), path('*report.txt'), path('*_kraken2_classified_reads.txt'), path(reads), path('*.unclassified{.,_}*'), emit: results

    

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def paired       = meta.single_end ? "" : "--paired"
    def classified   = meta.single_end ? "${prefix}.classified.fastq"   : "${prefix}.classified#.fastq"
    def unclassified = meta.single_end ? "${prefix}.unclassified.fastq" : "${prefix}.unclassified#.fastq"
    def classified_option = save_classified_output_fastqs ? "--classified-out ${classified}" : ""
    def unclassified_option = save_unclassified_output_fastqs ? "--unclassified-out ${unclassified}" : ""
    def readclassification_option = save_reads_assignment ? "--output ${prefix}_kraken2_classified_reads.txt" : "--output /dev/null"
    def compress_reads_command = (save_classified_output_fastqs || save_unclassified_output_fastqs) ? "pigz -p $task.cpus *.fastq" : ""

    """
    kraken2 \\
        --db $db \\
        --threads $task.cpus \\
        --report ${prefix}_kraken2_report.txt \\
        --gzip-compressed \\
        $unclassified_option \\
        $classified_option \\
        $readclassification_option \\
        $paired \\
        $args \\
        --report-minimizer-data \\
        --minimum-hit-groups 3 \\
        --confidence 0.05 \\
        --use-names \\
        $reads

    $compress_reads_command

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: \$(echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version //; s/ .*\$//')
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def paired       = meta.single_end ? "" : "--paired"
    def classified   = meta.single_end ? "${prefix}.classified.fastq.gz"   : "${prefix}.classified_1.fastq.gz ${prefix}.classified_2.fastq.gz"
    def unclassified = meta.single_end ? "${prefix}.unclassified.fastq.gz" : "${prefix}.unclassified_1.fastq.gz ${prefix}.unclassified_2.fastq.gz"
    def readclassification_option = save_reads_assignment ? "--output ${prefix}_kraken2_classified_reads.txt" : "--output /dev/null"
    def compress_reads_command = (save_classified_output_fastqs || save_unclassified_output_fastqs) ? "pigz -p $task.cpus *.fastq" : ""

    """
    touch ${prefix}.kraken2.report.txt
    if [ "$save_output_fastqs" == "true" ]; then
        touch $classified
        touch $unclassified
    fi
    if [ "$save_reads_assignment" == "true" ]; then
        touch ${prefix}_kraken2_classified_reads.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: \$(echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version //; s/ .*\$//')
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """

}