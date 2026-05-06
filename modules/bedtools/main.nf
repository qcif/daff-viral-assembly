process BEDTOOLS {
    tag "${sampleid}"
    label 'setting_3'
    publishDir { "${params.outdir}/${sampleid}/09_mapping_to_ref" }, mode: 'copy'

    input:
    tuple val(sampleid), path(ref), path(bam), path(bai), path(vcf_applied_fasta)

    output:
    path("${sampleid}_bcftools_masked_consensus.fasta")
    tuple val(sampleid), path("${sampleid}_bcftools_masked_consensus.fasta"), emit: bcftools_masked_consensus_fasta

    script:
    """
    bedtools genomecov -ibam ${bam} -bga > ${sampleid}_genomecov.bed
    awk '\$4==0 {print}' ${sampleid}_genomecov.bed > ${sampleid}_zero_coverage.bed
    bedtools maskfasta -fi ${vcf_applied_fasta} -bed ${sampleid}_zero_coverage.bed -fo ${sampleid}_bcftools_masked_consensus.fasta
    """
}