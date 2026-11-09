process COVSTATS {
    tag "$sampleid"
    label 'setting_1'

    publishDir {
        mode == 'reference'
            ? "${params.outdir}/${sampleid}/09_mapping_to_ref"
            : "${params.outdir}/${sampleid}/08_mapping_to_contigs"
    }, mode: 'copy'

    input:
    tuple val(sampleid),
          val(mode),
          path(bed),
          path(per_base),
          path(blast_results),
          path(bbsplit_stats),
          path(coverage),
          path(mapping_q)

    output:
    path("*_with_cov_stats.txt")
    tuple val(sampleid),
          path("*_with_cov_stats.txt"),
          emit: detections_summary

    script:
    def output_name = mode == 'reference' ? "${sampleid}_reference_with_cov_stats.txt" : "${sampleid}_contigs_with_cov_stats.txt"
    """
    if [ -s "${blast_results}" ] && [ -s "${bbsplit_stats}" ] && [ -s "${coverage}" ] && [ -s "${bed}" ] && [ -s "${mapping_q}" ] && [ -s "${per_base}" ]; then
        derive_coverage_stats.py \\
          --mode ${mode} \\
          --sample ${sampleid} \\
          --blastn_results ${blast_results} \\
          --bbsplit_stats ${bbsplit_stats} \\
          --coverage ${coverage} \\
          --bed ${bed} \\
          --mapping_quality ${mapping_q} \\
          --per_base ${per_base}
    else
        touch "${output_name}"
    fi
    """
}