process SAMTOOLS_MPILEUP {
    publishDir { "${params.outdir}/${sampleid}/08_mapping_to_contigs" }, mode: 'copy'
    tag "${sampleid}"
    label 'setting_10'

    input:
    tuple val(sampleid), path(ref), path(sam)

    output:
    tuple val(sampleid), path("${sampleid}_contig_pileup.txt"), path(ref), emit: pileup

    script:
    """
    #filter out unmapped reads -F 4
    samtools view -@ ${task.cpus} -Sb -F 4 ${sam} | samtools sort -@ ${task.cpus} -o ${sampleid}_contig_aln.sorted.bam
    samtools index ${sampleid}_contig_aln.sorted.bam
    samtools mpileup -aa -f ${ref} ${sampleid}_contig_aln.sorted.bam > ${sampleid}_contig_pileup.txt  
    """
}