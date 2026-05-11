/*
process FASTA2TABLE2 {
    tag "$sampleid"
    label 'setting_1'
    publishDir { "${params.outdir}/${sampleid}/09_mapping_to_ref" }, mode: 'copy'

    input:
    tuple val(sampleid), path(fasta), path(stats)
    
    output:
    file("${sampleid}_reference_with_cov_stats_final.txt") done
    tuple val(sampleid), file("${sampleid}_reference_with_cov_stats_final.txt"), emit: detections_summary_final

    script:
    """
    fasta2table.py --fasta ${fasta} --sample ${sampleid} --tophits ${stats} --mode reference
    """
}
*/
process FASTA2TABLE_REF {
    tag "$sampleid"
    label 'setting_1'
    publishDir { "${params.outdir}/${sampleid}/09_mapping_to_ref" }, mode: 'copy'

    input:
    tuple val(sampleid), path(stats), path(fasta), val(run_mode)

    output:
    // reference mode outputs
    path("${sampleid}_reference_with_cov_stats_final.txt")
    tuple val(sampleid), path("${sampleid}_reference_with_cov_stats_final.txt"), emit: detections_summary_final

    script:
    """
    fasta2table.py --fasta ${fasta} --sample ${sampleid} --tophits ${stats} --mode ${run_mode}
    """
}
