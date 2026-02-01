/*
process FASTP {
  tag "$sampleid"
  publishDir "${params.outdir}/${sampleid}/02_qtrimmed", mode: 'copy', pattern: '{*fastq.gz,*_fastp.log}'
  publishDir "${params.outdir}/${sampleid}/03_fastqc_trimmed", mode: 'copy', pattern: '{*html,*json}'
  label "setting_4"

  input:
    tuple val(sampleid), path(fastq1), path(fastq2)

  output:
    path("${sampleid}.fastp.html")
    path("${sampleid}.fastp.json")
    path("${sampleid}_fastp.log")
    tuple val(sampleid), path("${sampleid}_1_qtrimmed.fastq.gz"), path("${sampleid}_2_qtrimmed.fastq.gz"), emit: trimmed_fq
    path("${sampleid}.fastp.json"), emit: fastp_json

  script:
    """
    fastp \
    --in1 ${fastq1} \
    --in2 ${fastq2} \
    --out1 ${sampleid}_1_qtrimmed.fastq.gz \
    --out2 ${sampleid}_2_qtrimmed.fastq.gz \
    --cut_front \
    --cut_tail \
    --json ${sampleid}.fastp.json \
    --html ${sampleid}.fastp.html \
    --thread 6 \
    --detect_adapter_for_pe \
    --length_required 50 --average_qual 20
    2>&1 | tee ${sampleid}_fastp.log
    """
}
//Modified the nf-core module so it does not expect an adapter list. Might re-visit later to add that functionality back in.
//Also added --detect_adapter_for_pe --cut_front --cut_tail --length_required 50 --average_qual 20, might want to make these external arguments later.
*/
process FASTP {
    tag "$meta.id"
    label "setting_4"
    publishDir "${params.outdir}/$meta.id/02_qtrimmed", mode: 'copy', pattern: '{*fastq.gz}'
    publishDir "${params.outdir}/$meta.id/03_fastqc_trimmed", mode: 'copy', pattern: '{*fastp.html,*fastp.json}'

    conda "bioconda::fastp=0.23.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastp:0.23.4--h5f740d0_0' :
        'biocontainers/fastp:0.23.4--h5f740d0_0' }"

    input:
    tuple val(meta), path(reads)
    //path  adapter_fasta
    val   save_trimmed_fail
    val   save_merged

    output:
    path("*.fastp.fastq.gz")
    path("*.fastp.html")
    path("*.fastp.json")
    tuple val(meta), path('*.fastp.fastq.gz') , emit: reads
    tuple val(meta), path('*.json')           , emit: json
    tuple val(meta), path('*.html')           , emit: html
    tuple val(meta), path('*.log')            , emit: log
    path "versions.yml"                       , emit: versions
    tuple val(meta), path('*.fail.fastq.gz')  , optional:true, emit: reads_fail
    tuple val(meta), path('*.merged.fastq.gz'), optional:true, emit: reads_merged

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    //def adapter_list = adapter_fasta ? "--adapter_fasta ${adapter_fasta}" : ""
    def fail_fastq = save_trimmed_fail && meta.single_end ? "--failed_out ${prefix}.fail.fastq.gz" : save_trimmed_fail && !meta.single_end ? "--unpaired1 ${prefix}_1.fail.fastq.gz --unpaired2 ${prefix}_2.fail.fastq.gz" : ''
    // Added soft-links to original fastqs for consistent naming in MultiQC
    // Use single ended for interleaved. Add --interleaved_in in config.
    if ( task.ext.args?.contains('--interleaved_in') ) {
        """
        [ ! -f  ${prefix}.fastq.gz ] && ln -sf $reads ${prefix}.fastq.gz

        fastp \\
            --stdout \\
            --in1 ${prefix}.fastq.gz \\
            --thread $task.cpus \\
            --json ${prefix}.fastp.json \\
            --html ${prefix}.fastp.html \\
            --detect_adapter_for_pe \\
            --cut_front \\
            --cut_tail \\
            --length_required 50 --average_qual 20 \\
            $fail_fastq \\
            $args \\
            2> >(tee ${prefix}.fastp.log >&2) \\
        | gzip -c > ${prefix}.fastp.fastq.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
        END_VERSIONS
        """
    } else if (meta.single_end) {
        """
        [ ! -f  ${prefix}.fastq.gz ] && ln -sf $reads ${prefix}.fastq.gz

        fastp \\
            --in1 ${prefix}.fastq.gz \\
            --out1  ${prefix}.fastp.fastq.gz \\
            --thread $task.cpus \\
            --json ${prefix}.fastp.json \\
            --html ${prefix}.fastp.html \\
            $fail_fastq \\
            $args \\
            2> >(tee ${prefix}.fastp.log >&2)

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
        END_VERSIONS
        """
    } else {
        def merge_fastq = save_merged ? "-m --merged_out ${prefix}.merged.fastq.gz" : ''
        """
        [ ! -f  ${prefix}_1.fastq.gz ] && ln -sf ${reads[0]} ${prefix}_1.fastq.gz
        [ ! -f  ${prefix}_2.fastq.gz ] && ln -sf ${reads[1]} ${prefix}_2.fastq.gz
        fastp \\
            --in1 ${prefix}_1.fastq.gz \\
            --in2 ${prefix}_2.fastq.gz \\
            --out1 ${prefix}_1.fastp.fastq.gz \\
            --out2 ${prefix}_2.fastp.fastq.gz \\
            --json ${prefix}.fastp.json \\
            --html ${prefix}.fastp.html \\
            $fail_fastq \\
            $merge_fastq \\
            --thread $task.cpus \\
            --detect_adapter_for_pe \\
            $args \\
            2> >(tee ${prefix}.fastp.log >&2)

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
        END_VERSIONS
        """
    }
}