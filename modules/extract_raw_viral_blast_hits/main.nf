process EXTRACT_RAW_VIRAL_BLAST_HITS {
    tag "${sampleid}"
    label 'setting_7'
    //containerOptions "--bind ${file(params.taxdump)}"
    publishDir { "${params.outdir}/${sampleid}/07_annotation" }, mode: 'copy'

    input:
    tuple val(sampleid), path(blast_results), path(assembly_headers)
    val(mounted_db_dir)
    path(staged_db_dir)
    path(filter_terms)

    output:
    tuple val(sampleid), path("${sampleid}_megablast_top_viral_hits.txt"), emit: viral_blast_results
    path("${sampleid}_blastn.txt")


    script:
    def taxonkit_db = mounted_db_dir ?: staged_db_dir
    """
    cat ${blast_results} > ${sampleid}_blastn.txt
    filter_blast.py --blastn_results ${sampleid}_blastn.txt \
                    --sample_name ${sampleid} \
                    --taxonkit_database_dir ${taxonkit_db} \
                    --filter ${filter_terms} \
                    --assembly_headers ${assembly_headers}
    """
}