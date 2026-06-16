process SEQTK_SUBSEQ {
    tag "${sampleid}"
    label 'setting_1'

    input:
    //tuple val(sampleid), path(viral_contig_ids), path(other_contig_ids), path(contigs)
    tuple val(sampleid), path(blast_results), path(fasta)
    output:
    tuple val(sampleid), path("${sampleid}_candidate_viral_contigs.fasta"), emit: viral_candidate_fasta
    tuple val(sampleid), path("${sampleid}_other_contigs.fasta"), emit: other_fasta

    script:
    """
    cut -f2 ${blast_results} | sed '1d' | sed 's/ //g' | sort | uniq > ${sampleid}_viral_candidate_contig_ids.txt
    grep '>' ${fasta} | sed 's/>//g' | sed 's/ /_/g' | sort | uniq > ${sampleid}_all_contigs.txt
    grep -F -x -v -f ${sampleid}_viral_candidate_contig_ids.txt \
      ${sampleid}_all_contigs.txt \
      > ${sampleid}_other_contig_ids.txt || touch ${sampleid}_other_contig_ids.txt
    if [[ ! -s ${sampleid}_viral_candidate_contig_ids.txt ]]; then
      touch ${sampleid}_candidate_viral_contigs.fasta
    else
      seqtk subseq ${fasta} ${sampleid}_viral_candidate_contig_ids.txt > ${sampleid}_candidate_viral_contigs.fasta
      seqtk subseq ${fasta} ${sampleid}_other_contig_ids.txt > ${sampleid}_other_contigs.fasta
    fi  
    seqtk subseq ${fasta} ${sampleid}_viral_candidate_contig_ids.txt > ${sampleid}_candidate_viral_contigs.fasta
    seqtk subseq ${fasta} ${sampleid}_other_contig_ids.txt > ${sampleid}_other_contigs.fasta
    """
}