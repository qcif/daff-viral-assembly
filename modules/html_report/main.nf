process HTML_REPORT {
    publishDir { "${params.outdir}/${sampleid}/11_report" }, mode: 'copy', overwrite: true
    label 'setting_4'

    input:
    tuple val(sampleid), path(raw_fastqc), path(filtered_fastqc), path(fastp), path(fasta), path(summary_known_viruses), path(kaiju_summary), path(kraken_summary), path(detections_summary), path(ref_mapping_summary), path(consensus), path(bam), path(bai), path(novel_virus_summary), path(novel_support_summary), path(blast_contig2ref),path(orfs), path(hmmscan), path(diamond_summary),
    path(timestamp),
    path(qcreport_html),
    path(qcreport_txt),
    path(configyaml),
    path(samplesheet)

    output:
    path("*"), optional: true
    path(raw_fastqc)
    path(filtered_fastqc)
    path(qcreport_html)
    path(bam)
    path(bai)

    script:
    analyst_name = params.analyst_name ? params.analyst_name.replaceAll(/ /, '_') : "unknown"
    facility = params.facility ? params.facility.replaceAll(/ /, '_') : "unknown"
    """
    #cp ${qcreport_html} .
    cp ${params.tool_versions} versions.yml
    cp ${params.default_params} default_params.yml
    cp ${params.filter_terms} filterKeyWords.txt

    build_report.py --samplesheet ${samplesheet} --result_dir . --params_file ${configyaml} --analyst ${analyst_name} --facility ${facility} --versions versions.yml --default_params_file default_params.yml
    """
}
