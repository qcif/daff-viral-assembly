process SEQTK_SUBSEQ {
    tag "${sampleid}"
    label "setting_1"

    input:
    tuple val(sampleid), path(viral_contig_ids), path(other_contig_ids), path(contigs)
    
    output:
    tuple val(sampleid), path("${sampleid}_candidate_viral_contigs.fasta"), emit: viral_candidate_fasta
    tuple val(sampleid), path("${sampleid}_other_contigs.fasta"), emit: other_fasta

    script:
    """
    if [[ ! -s ${viral_contig_ids} ]]; then
      touch ${sampleid}_candidate_viral_contigs.fasta
    else
      seqtk subseq ${contigs} ${viral_contig_ids} > ${sampleid}_candidate_viral_contigs.fasta
      seqtk subseq ${contigs} ${other_contig_ids} > ${sampleid}_other_contigs.fasta
    fi
    """
}