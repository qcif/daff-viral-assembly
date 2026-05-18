process SUMMARISE_RESULTS {
    tag "${sampleid}"
    label 'setting_3'
    publishDir { "${params.outdir}/${sampleid}/10_results_summary" }, mode: 'copy', pattern: '{*summary_viral_results.tsv}'
    publishDir { "${params.outdir}/${sampleid}/10_results_summary" }, mode: 'copy', pattern: '{*novel_virus_candidates.tsv}'
    publishDir { "${params.outdir}/${sampleid}/10_results_summary" }, mode: 'copy', pattern: '{*novel_evidence_summary.txt}'
    publishDir { "${params.outdir}/${sampleid}/07_annotation" }, mode: 'copy', pattern: '{*hmm_domain_summary_counts.tsv}'

    input:
    tuple val(sampleid), path(kraken_results), path(kaiju_results), path(blast), path(hmmscan), path(map2ref), path(contigs), path(genomad), path(blast_novel), path(diamond_results), path(taxonomy)

    output:
    path("${sampleid}_summary_viral_results.tsv")
    path("${sampleid}_hmm_domain_summary_counts.tsv")
    path("${sampleid}_novel_virus_candidates.tsv")
    path("${sampleid}_novel_evidence_summary.txt")
    tuple val(sampleid), path("${sampleid}_hmm_domain_summary_counts.tsv"), emit: domain_count
    tuple val(sampleid), path("${sampleid}_summary_viral_results.tsv"), emit: summary_known_viruses
    tuple val(sampleid), path("${sampleid}_novel_virus_candidates.tsv"), emit: novel_virus_candidates
    tuple val(sampleid), path("${sampleid}_novel_evidence_summary.txt"), emit: novel_support_summary

    script:
    """
    viral_results_summary.py \\
      --sample_name ${sampleid} \\
      --kaiju ${kaiju_results} \\
      --kraken ${kraken_results} \\
      --blast ${blast} \\
      --hmmscan ${hmmscan} \\
      --map2ref ${map2ref} \\
      --genomad ${genomad} \\
      --blast_novel ${blast_novel} \\
      --diamond ${diamond_results} \\
      --taxonomy ${taxonomy} \\
      --min-reads 2000 \\
      --fasta ${contigs}
    """
}