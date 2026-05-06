process HMMSCAN {
    tag "${sampleid}"
    label "setting_20"
    publishDir { "${params.outdir}/${sampleid}/07_annotation" }, mode: 'copy'
    containerOptions "${params.bindOptions}"

    input:
    tuple val(sampleid), path(fasta)
    val(hmmer_db)
    
    output:
    file "${sampleid}_orfs.fasta"
    file "${sampleid}_hmmscan*_output.txt"
    tuple val(sampleid), path("${sampleid}_hmmscan_per_target_output.txt"), emit: hmmscan_preds
    tuple val(sampleid), path("${sampleid}_hmmscan_per_domain_output.txt"), emit: hmmscan_domain_preds

    script:
    """
    hmmscan --cpu ${task.cpus} \\
            --domtblout ${sampleid}_hmmscan_per_domain_output.txt \\
            --tblout ${sampleid}_hmmscan_per_target_output.txt \\
            --pfamtblout ${sampleid}_hmmscan_succinct_output.txt \\
            ${hmmer_db} ${fasta} \\
            > ${sampleid}_hmmscan.log 2>&1
    """
}