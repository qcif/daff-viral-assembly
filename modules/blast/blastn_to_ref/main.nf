process BLAST_BLASTN_TO_REF {
    tag "${sampleid}"
    label 'setting_20'

    input:
    tuple val(sampleid), path(assembly), path(reference)

    output:
    tuple val(sampleid), path("${sampleid}*_blastn.bls"), emit: blast_results

    script:
    """
    makeblastdb -in ${reference} -parse_seqids -dbtype nucl
    blastn -query ${assembly} \\
      -db ${reference} \\
      -out ${sampleid}_contig_vs_refs_blastn.bls \\
      -evalue 1e-3 \\
      -num_threads 2 \\
      -max_target_seqs 1
    """
}