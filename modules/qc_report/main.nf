process QC_REPORT {
    label 'setting_1'
    publishDir "${params.outdir}/02_qc_report", mode: 'copy', overwrite: true

    input:
    path multiqc_files

    output:
    path("run_qc_report_*txt")
    path("run_qc_report_*html")
    path("run_qc_report_*html"), emit: qc_report_html
    path("run_qc_report_*txt"), emit: qc_report_txt

    script:
    """
    seq_run_qc_report.py --qfiltered_reads_threshold ${params.qfiltered_reads_threshold} --clean_reads_threshold ${params.clean_reads_threshold}
    """
}