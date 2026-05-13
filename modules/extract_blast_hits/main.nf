process EXTRACT_BLAST_HITS {
    tag "${sampleid}"
    label 'setting_1'

    input:
    tuple val(sampleid), path(blast_results), val(target_organism), val(target_gene), val(target_size)

    output:
    tuple val(sampleid), path("${sampleid}*_megablast_top_hits_tmp.txt"), emit: topblast
    tuple val(sampleid), path("${sampleid}_reference_match.fasta"), emit: reference_fasta_files
    tuple val(sampleid), path("${sampleid}_final_polished_consensus_match.fasta"), emit: consensus_fasta_files

    script:
    target_organism_str = (target_organism instanceof List)
    ? "\"${target_organism.join('|')}\""
    : "\"${target_organism}\""
    """
    if [[ \$(wc -l < *_megablast_top_10_hits.txt) -ge 2 ]]
      then
        select_top_blast_hit.py --sample_name ${sampleid} --blastn_results ${sampleid}*_top_10_hits.txt --target_organism ${target_organism_str} --taxonkit_database_dir ${params.taxdump}

        # extract segment of consensus sequence that align to reference
        awk  -F  '\\t' 'NR>1 { printf ">%s\\n%s\\n",\$2,\$23 }' ${sampleid}*_top_hits_tmp.txt | sed 's/-//g' > ${sampleid}_final_polished_consensus_match.fasta

        # extract segment of reference that align to consensus sequence
        awk  -F  '\\t' 'NR>1 { printf ">%s_%s\\n%s\\n",\$2,\$4,\$24 }' ${sampleid}*_top_hits_tmp.txt | sed 's/-//g' > ${sampleid}_reference_match.fasta
    else
        echo "No hits found for ${sampleid} in the blast results. Skipping the extraction of consensus and reference fasta files." >&2
        touch ${sampleid}_final_polished_consensus_match.fasta
        touch ${sampleid}_reference_match.fasta
        touch ${sampleid}_megablast_top_hits_tmp.txt
    fi
    """
}