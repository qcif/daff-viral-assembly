/*
process FASTA2TABLE {
    tag "$sampleid"
    label 'setting_1'

    input:
    tuple val(sampleid), path(tophits), path(fasta)
    output:
    file("${sampleid}_megablast_top_viral_hits_with_contigs.txt") done
    file("${sampleid}_megablast_top_viral_hits_filtered_with_contigs.txt") done
    tuple val(sampleid), file("${sampleid}_ref_ids_to_retrieve.txt"), emit: ref_ids
    tuple val(sampleid), file("${sampleid}*_filtered_viral_contigs.fasta"), emit: contig_fasta
    tuple val(sampleid), file("${sampleid}_megablast_top_viral_hits_filtered_with_contigs.txt"), emit: blast_results
    tuple val(sampleid), file("${sampleid}_megablast_top_viral_hits_with_contigs.txt"), emit: blast_results2

    script:
    """
    fasta2table.py --fasta ${fasta} --sample ${sampleid} --tophits ${tophits} --mode contigs
    cut -f2 ${sampleid}_megablast_top_viral_hits_filtered_with_contigs.txt | sed '1d' | sed 's/ //g' | sort | uniq > ${sampleid}_contig_ids_to_retrieve.txt
    """
}

*/
process FASTA2TABLE_CONTIGS {
    tag "$sampleid"
    label 'setting_1'

    input:
    tuple val(sampleid), path(stats), path(fasta), val(run_mode)

    output:
    // contigs mode outputs
    tuple val(sampleid), path("${sampleid}_ref_ids_to_retrieve.txt"), emit: ref_ids
    tuple val(sampleid), path("${sampleid}*_filtered_viral_contigs.fasta"), emit: contig_fasta
    tuple val(sampleid), path("${sampleid}_megablast_top_viral_hits_filtered_with_contigs.txt"), emit: blast_results
    tuple val(sampleid), path("${sampleid}_megablast_top_viral_hits_with_contigs.txt"), emit: blast_results2

    script:
    """
    fasta2table.py --fasta ${fasta} --sample ${sampleid} --tophits ${stats} --mode ${run_mode}
    """
}
