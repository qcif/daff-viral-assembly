
process IDENTIFY_ERRORS {
    tag "${sampleid}"
    label 'setting_2'
    publishDir { "${params.outdir}/${sampleid}/08_mapping_to_contigs" }, mode: 'copy'

    input:
    tuple val(sampleid), path(pileup), path(ref)

    output:
    tuple val(sampleid), path("${sampleid}_trim_coords.tsv"), path(ref), emit: trimmed_coords

    script:
    """
    trim_contig_ends.py ${pileup} > ${sampleid}_trim_coords.tsv
    """
}