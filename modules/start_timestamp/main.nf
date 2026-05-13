process START_TIMESTAMP {
    publishDir "${params.outdir}/01_pipeline_logs", mode: 'copy', overwrite: true
    cache false

    output:
    path "*nextflow_start_timestamp.txt"
    path("*nextflow_start_timestamp.txt"), emit: timestamp

    script:
    """
    START_TIMESTAMP=\$(date "+%Y%m%d%H%M%S")
    echo "\$START_TIMESTAMP" > "\${START_TIMESTAMP}_nextflow_start_timestamp.txt"
    """
}