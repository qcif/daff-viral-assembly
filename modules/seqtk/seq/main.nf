process SEQTK_SEQ {
  tag "${sampleid}"
  label "setting_1"
  publishDir { "${params.outdir}/${sampleid}/06_assembly" }, mode: 'copy', pattern: '{*_spades_scaffolds.fasta}'

  input:
    tuple val(sampleid), path(assembly)
  output:
    path("${sampleid}_spades_scaffolds.fasta")
    tuple val(sampleid), path("${sampleid}_spades_scaffolds.fasta"), emit: filt_fasta
    tuple val(sampleid), path("${sampleid}_spades_scaffolds_headers.txt"), emit: filt_headers

  script:
  """
  seqtk seq -L 150 ${assembly} > ${sampleid}_scaffolds_filt.fasta

  grep '>' ${sampleid}_scaffolds_filt.fasta | sed 's/>NODE/CONTIG/g' > ${sampleid}_spades_scaffolds_headers.txt
  sed -E 's/^>NODE_([0-9]+).*/>CONTIG_\\1/' ${sampleid}_scaffolds_filt.fasta > ${sampleid}_spades_scaffolds.fasta
  """
}