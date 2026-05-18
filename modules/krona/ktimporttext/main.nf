process KRONA_KTIMPORTTEXT {
    publishDir { "${params.outdir}/${sampleid}/05_read_classification" }, mode: 'link'
    label 'setting_2'
    tag "${sampleid}"

    input:
    tuple val(sampleid), path(krona_input)

    output:
    file "${sampleid}_kaiju_krona.html"

    script:
    """
    ktImportText -o ${sampleid}_kaiju_krona.html ${krona_input}
    """
}