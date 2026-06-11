process SAMTOOLS_REF {
    publishDir { "${params.outdir}/${sampleid}/09_mapping_to_ref" }, mode: 'copy'
    tag "${sampleid}"
    label 'setting_10'

    input:
    tuple val(sampleid), path(ref), path(sam)

    output:
    path "${sampleid}_ref_aln.sorted.bam"
    path "${sampleid}_ref_aln.sorted.bam.bai"
    path "${sampleid}_samtools_consensus_from_ref.fasta"
    path "${sampleid}_ref_coverage.txt"
    tuple val(sampleid), path(ref), path("${sampleid}_ref_aln.sorted.bam"), path("${sampleid}_ref_aln.sorted.bam.bai"), emit: sorted_bam
    tuple val(sampleid), path("${sampleid}_ref_coverage.txt"), emit: coverage
    tuple val(sampleid), path("${sampleid}_ref_mapq.txt"), emit: mapping_quality

    script:
    """
    #filter out unmapped reads -F 4
    samtools view -@ ${task.cpus} -Sb -F 4 ${sam} | samtools sort -@ ${task.cpus} -o ${sampleid}_ref_aln.sorted.bam
    samtools index ${sampleid}_ref_aln.sorted.bam
    samtools coverage ${sampleid}_ref_aln.sorted.bam  > ${sampleid}_ref_coverage.txt
    samtools view -@ ${task.cpus} ${sampleid}_ref_aln.sorted.bam | awk '{mapq[\$3]+=\$5; count[\$3]++} END {for (chr in mapq) printf "%s\\t%.2f\\n", chr, mapq[chr]/count[chr]}' > ${sampleid}_ref_mapq.txt
    samtools consensus -@ ${task.cpus} -f fastq -a -A ${sampleid}_ref_aln.sorted.bam -o ${sampleid}_consensus.fastq 
    samtools consensus -@ ${task.cpus} -f fasta -a -A ${sampleid}_ref_aln.sorted.bam -o ${sampleid}_samtools_consensus_from_ref.fasta
    """
}
