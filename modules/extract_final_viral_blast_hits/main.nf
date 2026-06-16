process EXTRACT_FINAL_VIRAL_BLAST_HITS {
    tag "${sampleid}"
    label 'setting_7'
    publishDir { "${params.outdir}/${sampleid}/07_annotation" }, mode: 'copy'

    input:
    tuple val(sampleid), path(blast_results), path(assembly_headers)
    path(taxonkit_db)

    output:
    tuple val(sampleid), path("${sampleid}_megablast_top_viral_hits.txt"), emit: viral_blast_results
    path("${sampleid}_blastn.txt")

    script:
    """
    cat ${blast_results} > ${sampleid}_blastn.txt
    filter_blast.py --blastn_results ${sampleid}_blastn.txt \
                    --sample_name ${sampleid} \
                    --taxonkit_database_dir ${taxonkit_db} \
                    --filter ${params.filter_terms} \
                    --assembly_headers ${assembly_headers}
    """
}