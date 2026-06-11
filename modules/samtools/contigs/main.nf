process SAMTOOLS_CONTIGS {
    publishDir { "${params.outdir}/${sampleid}/08_mapping_to_contigs" }, mode: 'copy'
    tag "${sampleid}"
    label 'setting_10'

    input:
    tuple val(sampleid), path(ref), path(sam)

    output:
    path "${sampleid}_contig_aln.sorted.bam"
    path "${sampleid}_contig_aln.sorted.bam.bai"
    path "${sampleid}_contig_coverage.txt"
    tuple val(sampleid), path(ref), path("${sampleid}_contig_aln.sorted.bam"), path("${sampleid}_contig_aln.sorted.bam.bai"), emit: sorted_bam
    tuple val(sampleid), path("${sampleid}_contig_coverage.txt"), emit: coverage
    tuple val(sampleid), path("${sampleid}_contig_mapq.txt"), emit: mapping_quality

    script:
    """
    #filter out unmapped reads -F 4
    samtools view -@ ${task.cpus} -Sb -F 4 ${sam} | samtools sort -@ ${task.cpus} -o ${sampleid}_contig_aln.sorted.bam
    samtools index ${sampleid}_contig_aln.sorted.bam
    samtools coverage ${sampleid}_contig_aln.sorted.bam  > ${sampleid}_contig_coverage.txt
    samtools view -@ ${task.cpus} ${sampleid}_contig_aln.sorted.bam | awk '{mapq[\$3]+=\$5; count[\$3]++} END {for (chr in mapq) printf "%s\\t%.2f\\n", chr, mapq[chr]/count[chr]}' > ${sampleid}_contig_mapq.txt
    """
}