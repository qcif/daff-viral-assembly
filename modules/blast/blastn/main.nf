process BLAST_BLASTN {
    tag "${sampleid}"
    label 'setting_8'
    //containerOptions "--bind ${file(params.blastn_db).parent}"

    input:
    tuple val(sampleid), path(assembly)
    tuple path(db_dir), val(db_name)
    
    output:
    tuple val(sampleid), path("${sampleid}*_blastn.bls"), emit: blast_results

    script:
    //def blastdb_dir  = file(db).parent
    //def blastdb_name = file(db).name
    def blastoutput = assembly.getBaseName() + "_blastn.bls"
    """
    export BLASTDB=${db_dir}
    blastn -query ${assembly} \\
        -db ${db_dir}/${db_name} \\
        -out ${blastoutput} \\
        -evalue 1e-3 \\
        -num_threads ${params.blast_threads} \\
        -outfmt '6 qseqid sgi sacc length pident mismatch gapopen qstart qend qlen sstart send slen sstrand evalue bitscore qcovhsp stitle staxids qseq sseq sseqid qcovs qframe sframe' \\
        -max_target_seqs 5
    """
}