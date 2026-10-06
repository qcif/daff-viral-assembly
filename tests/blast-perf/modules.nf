// Processes for the batch-vs-per-query MEGABLAST performance test.
// Shared by batch.nf, per_query.nf and node_up.nf.

params.blast_container = "quay.io/biocontainers/blast:2.16.0--h66d330f_4"

process NODE_UP {
    tag "node_up"
    container params.blast_container
    publishDir params.outdir, mode: 'copy'
    cache false

    input:
    val db_path

    output:
    val true, emit: done
    path "node_up.node.txt", emit: node

    script:
    """
    hostname > node_up.node.txt
    nproc >> node_up.node.txt
    free -g >> node_up.node.txt
    blastdbcmd -db ${db_path} -info >> node_up.node.txt
    """
}

process BLAST_BATCH {
    tag "batch"
    container params.blast_container
    publishDir params.outdir, mode: 'copy'
    cache false

    input:
    path queries
    val db_path
    val ready

    output:
    path "batch.bls", emit: bls
    path "batch.timing", emit: timing
    path "batch.node.txt", emit: node

    script:
    """
    hostname > batch.node.txt
    nproc >> batch.node.txt
    free -g >> batch.node.txt

    start=\$(date +%s)
    blastn -task megablast -query ${queries} -db ${db_path} \\
        -evalue 1e-3 -max_target_seqs 5 -num_threads ${task.cpus} \\
        -outfmt '6 qseqid sgi sacc length pident mismatch gapopen qstart qend qlen sstart send slen sstrand evalue bitscore qcovhsp stitle staxids qseq sseq sseqid qcovs qframe sframe' \\
        -out batch.bls
    end=\$(date +%s)

    echo "start=\${start}" > batch.timing
    echo "end=\${end}" >> batch.timing
    echo "seconds=\$(awk -v s=\${start} -v e=\${end} 'BEGIN { print e - s }')" >> batch.timing
    """
}

process BLAST_SINGLE {
    tag "${query_id}"
    container params.blast_container
    publishDir params.outdir, mode: 'copy'
    cache false

    input:
    tuple val(query_id), path(query_fasta)
    val db_path
    val ready

    output:
    path "${query_id}.bls", emit: bls
    path "${query_id}.timing", emit: timing
    path "${query_id}.node.txt", emit: node

    script:
    """
    hostname > ${query_id}.node.txt
    nproc >> ${query_id}.node.txt
    free -g >> ${query_id}.node.txt

    start=\$(date +%s)
    blastn -task megablast -query ${query_fasta} -db ${db_path} \\
        -evalue 1e-3 -max_target_seqs 5 -num_threads ${task.cpus} \\
        -outfmt '6 qseqid sgi sacc length pident mismatch gapopen qstart qend qlen sstart send slen sstrand evalue bitscore qcovhsp stitle staxids qseq sseq sseqid qcovs qframe sframe' \\
        -out ${query_id}.bls
    end=\$(date +%s)

    echo "start=\${start}" > ${query_id}.timing
    echo "end=\${end}" >> ${query_id}.timing
    echo "seconds=\$(awk -v s=\${start} -v e=\${end} 'BEGIN { print e - s }')" >> ${query_id}.timing
    """
}
