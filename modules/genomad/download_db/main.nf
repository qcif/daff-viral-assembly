process GENOMAD_DOWNLOAD_DB {
    publishDir "${params.databases}", mode: 'copy'
    label "setting_2"

    output:
    path "databases/genomad_db", emit: db

    script:
    """
    mkdir -p databases
    genomad download-database databases
    """
}
