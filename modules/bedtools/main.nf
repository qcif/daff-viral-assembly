process BEDTOOLS {
    tag "${sampleid}"
    label 'setting_7'
    publishDir { "${params.outdir}/${sampleid}/09_mapping_to_ref" }, mode: 'copy'

    input:
    tuple val(sampleid), path(ref), path(bam), path(bai)

    output:
    tuple val(sampleid), path(ref), path(bam), path(bai), path("${sampleid}_zero_coverage.bed"), emit: bed_results
    script:
    """
    bedtools genomecov -ibam ${bam} -bga > ${sampleid}_genomecov.bed
    awk '\$4==0 {print}' ${sampleid}_genomecov.bed > ${sampleid}_zero_coverage.bed
    """
}