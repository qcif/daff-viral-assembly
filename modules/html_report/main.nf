process HTML_REPORT {
    publishDir { "${params.outdir}/${sampleid}/11_report" }, mode: 'copy', overwrite: true
    label 'setting_4'

    input:
    tuple val(sampleid), path(raw_fastqc), path(filtered_fastqc), path(fastp), path(fasta), path(summary_known_viruses), path(kaiju_summary), path(kraken_summary), path(detections_summary), path(ref_mapping_summary), path(consensus), path(bam), path(bai), path(novel_support_summary), path(blast_contig2ref),path(orfs), path(hmmscan), path(diamond_summary), path(novel_contig_summary),
    path(timestamp),
    path(qcreport_html),
    path(qcreport_txt),
    path(configyaml),
    path(versions_yml),
    path(default_params_yml),
    path(filter_terms_txt),
    path(samplesheet)

    output:
    path("${sampleid}_report.html")
    path("report_context.json")
    path(raw_fastqc)
    path(filtered_fastqc)
    path(qcreport_html)
    path(bam)
    path(bai)

    script:
    analyst_name = params.analyst_name ? params.analyst_name.replaceAll(/ /, '_') : "unknown"
    facility = params.facility ? params.facility.replaceAll(/ /, '_') : "unknown"
    """
    set +e
    build_report.py --samplesheet ${samplesheet} --result_dir . --params_file ${configyaml} --analyst ${analyst_name} --facility ${facility} --versions versions.yml --default_params_file default_params.yml
    
    status=\$?

    echo "BUILD_REPORT_EXIT_STATUS=\$status"
    echo "FILES AFTER REPORT:"
    ls -lah

    exit \$status
    """
}
