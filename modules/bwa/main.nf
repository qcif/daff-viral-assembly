process BWA {
  tag "${sampleid}"
  label "setting_21"

  input:
    tuple val(sampleid), path(contigs), path(fastq1), path(fastq2)

  output:
    tuple val(sampleid), path(contigs), file("${sampleid}_contig_aln.sam"), emit: contig_aligned_sam

  script:
  """
  bwa index ${contigs}
  bwa mem -t 4 ${contigs} $fastq1 $fastq2 > ${sampleid}_contig_aln.sam 2>> ${sampleid}_mapping.log
  """
}


/*
process REALIGN {
  tag "${sampleid}"
  label "setting_21"

  input:
    tuple val(sampleid), path(contigs), path(fastq1), path(fastq2)

  output:
    tuple val(sampleid), path(contigs), file("${sampleid}_contig_realn.sam"), emit: contig_aligned_sam

  script:
  """
  bwa index ${contigs}
  bwa mem -t 2 ${contigs} $fastq1 $fastq2 > ${sampleid}_contig_realn.sam 2>> ${sampleid}_mapping.log
  """
}
*/