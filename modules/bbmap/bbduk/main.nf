// original process stated explicitly inputs and outputs.
/*
process BBDUK { 
  tag "${sampleid}"
  label "setting_22"
  containerOptions "${bindOptions}"
  publishDir "${params.outdir}/${sampleid}/04_cleaned", mode: 'copy'

  input:
    tuple val(sampleid), path(fastq1), path(fastq2), path(ref)
  output:
    path("${sampleid}_rRNA_reads.log")
    path("${sampleid}_non_rRNA_fwd.fastq.gz")
    path("${sampleid}_non_rRNA_rev.fastq.gz")
    path("${sampleid}_rRNA_reads.log"), emit: bbduk_stats
    tuple val(sampleid), path("${sampleid}_non_rRNA_fwd.fastq.gz"), path("${sampleid}_non_rRNA_rev.fastq.gz"), emit: bbduk_filtered_fq

  script:
  """
  bbduk.sh -Xmx10g in=${fastq1} \
                   in2=${fastq2} \
                   out=${sampleid}_non_rRNA_fwd.fastq.gz \
                   out2=${sampleid}_non_rRNA_rev.fastq.gz \
                   outm=${sampleid}_rRNA_fwd.fastq.gz \
                   outm2=${sampleid}_rRNA_rev.fastq.gz \
                   k=31 ref=${ref} \
                   2>${sampleid}_rRNA_reads.log
  """
}
*/

//might want to specify parameter k=31 outside of process in the future
process BBMAP_BBDUK {
    tag "$meta.id"
    label 'setting_22'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5a/5aae5977ff9de3e01ff962dc495bfa23f4304c676446b5fdf2de5c7edfa2dc4e/data' :
        'community.wave.seqera.io/library/bbmap_pigz:07416fe99b090fa9' }"
    publishDir "${params.outdir}/${meta.id}/04_cleaned", mode: 'copy'

    input:
    tuple val(meta), path(reads)
    path(contaminants)

    output:
    //path("${meta.id}_non_rRNA_1.fastq.gz")
    //path("${meta.id}_non_rRNA_2.fastq.gz")
    path("${meta.id}_bbduk.log")
    tuple val(meta), path('*.fastq.gz'), emit: reads
    tuple val(meta), path('*.log')     , emit: log
    path('*.log')                      , emit: log2
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def raw      = meta.single_end ? "in=${reads[0]}" : "in1=${reads[0]} in2=${reads[1]}"
    def trimmed  = meta.single_end ? "out=${prefix}.fastq.gz" : "out1=${prefix}_non_rRNA_1.fastq.gz out2=${prefix}_non_rRNA_2.fastq.gz"
    def contaminants_fa = contaminants ? "ref=$contaminants" : ''
    """
    bbduk.sh \\
        -Xmx${task.memory.toGiga()}g \\
        $raw \\
        $trimmed \\
        k=31 \\
        threads=${task.cpus} \\
        $args \\
        $contaminants_fa \\
        &>${prefix}_bbduk.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh | grep -v "]")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output_command  = meta.single_end ? "echo '' | gzip > ${prefix}.fastq.gz" : "echo '' | gzip > ${prefix}_1.fastq.gz ; echo '' | gzip > ${prefix}_2.fastq.gz"
    """
    touch ${prefix}.bbduk.log
    $output_command

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh | grep -v "]")
    END_VERSIONS
    """
}