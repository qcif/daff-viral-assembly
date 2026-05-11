process TRIM_ENDS {
    tag "${sampleid}"
    label 'setting_2'
    publishDir { "${params.outdir}/${sampleid}/08_mapping_to_contigs" }, mode: 'copy'

    input:
    tuple val(sampleid), path(trimmed_coords), path(ref)

    output:
    tuple val(sampleid), path("${sampleid}_contigs.trimmed.fa"), emit: trimmed_contigs

    script:
    """
    samtools faidx ${ref}
    while read ctg L R; do samtools faidx ${ref} "\${ctg}:\${L}-\${R}" | sed "1s/.*/>\${ctg}/" >> ${sampleid}_contigs.trimmed.fa; done < ${trimmed_coords}
    """
}