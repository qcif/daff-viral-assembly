process SUMMARISE_READ_CLASSIFICATION {
    tag "${sampleid}"
    label 'setting_1'
    publishDir { "${params.outdir}/${sampleid}/05_read_classification" }, mode: 'copy'

    input:
    tuple val(sampleid), path(kaiju_results), path(kraken2_results), path(stats)
    val(mounted_db_dir)
    path(staged_db_dir)
    path(filter_terms)

    output:
    path("${sampleid}_kaiju_summary.txt")
    path("${sampleid}_kraken_summary.txt")
    tuple val(sampleid), path("${sampleid}_kaiju_summary.txt"), emit: kaiju_summary
    tuple val(sampleid), path("${sampleid}_kraken_summary.txt"), emit: kraken_summary

    script:
    def taxonkit_db = mounted_db_dir ?: staged_db_dir
    """
    filter_classification_results.py --kaiju ${kaiju_results} --sample_name ${sampleid} --kraken2 ${kraken2_results} --taxonkit_database_dir ${taxonkit_db} --stats ${stats} --filter ${filter_terms}
    """
}