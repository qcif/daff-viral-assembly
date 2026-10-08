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

process BLAST_PASS {
    tag "${pass}_${threads}"
    container params.blast_container
    publishDir params.outdir, mode: 'copy'
    cache false

    input:
    path queries
    val db_path
    val threads
    val pass
    val ready

    output:
    val true, emit: done
    path "${pass}_${threads}.bls", emit: bls
    path "${pass}_${threads}.timing", emit: timing
    path "${pass}_${threads}.sample", emit: sample
    path "${pass}_${threads}.node.txt", emit: node

    script:
    """
    hostname > ${pass}_${threads}.node.txt
    nproc >> ${pass}_${threads}.node.txt
    free -g >> ${pass}_${threads}.node.txt
    cat /sys/fs/cgroup/cpu.max >> ${pass}_${threads}.node.txt 2>/dev/null \\
        || echo "cpu.max=n/a" >> ${pass}_${threads}.node.txt
    cat /sys/fs/cgroup/memory.max >> ${pass}_${threads}.node.txt 2>/dev/null \\
        || echo "memory.max=n/a" >> ${pass}_${threads}.node.txt

    start=\$(date +%s)
    blastn -task megablast -query ${queries} -db ${db_path} \\
        -evalue 1e-3 -max_target_seqs 5 -num_threads ${threads} \\
        -outfmt '6 qseqid sgi sacc length pident mismatch gapopen qstart qend qlen sstart send slen sstrand evalue bitscore qcovhsp stitle staxids qseq sseq sseqid qcovs qframe sframe' \\
        -out ${pass}_${threads}.bls &
    bpid=\$!

    # Sample blastn's cumulative CPU time (utime+stime, field 14+15 of
    # /proc/<pid>/stat, in clock ticks) and cumulative disk sectors read
    # (field 6 of /proc/diskstats for nvme0n1) once a second. Deltas are
    # computed by summarise_cores.py, not here. Neither file exists in the
    # local smoke test (no real blastn child stat quirks, no NVMe device),
    # so each read falls back to 0 rather than failing the task.
    while kill -0 \$bpid 2>/dev/null; do
        t=\$(date +%s)
        cpu=\$(awk '{print \$14+\$15}' /proc/\$bpid/stat 2>/dev/null || echo 0)
        sectors=\$(awk '\$3=="nvme0n1"{print \$6}' /proc/diskstats 2>/dev/null)
        echo "\${t} \${cpu:-0} \${sectors:-0}" >> ${pass}_${threads}.sample
        sleep 1
    done
    wait \$bpid
    end=\$(date +%s)

    echo "start=\${start}" > ${pass}_${threads}.timing
    echo "end=\${end}" >> ${pass}_${threads}.timing
    echo "threads=${threads}" >> ${pass}_${threads}.timing
    echo "seconds=\$(awk -v s=\${start} -v e=\${end} 'BEGIN { print e - s }')" >> ${pass}_${threads}.timing
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

// --- scale test (scale-test.md): ~30k queries, batch vs. split-into-10 ---

process SCALE_BATCH {
    tag "scale_batch"
    container params.blast_container
    publishDir params.outdir, mode: 'copy'
    cache false

    input:
    path queries
    val db_path
    val ready

    output:
    path "scale_batch.bls", emit: bls
    path "scale_batch.timing", emit: timing
    path "scale_batch.node.txt", emit: node

    script:
    """
    hostname > scale_batch.node.txt
    nproc >> scale_batch.node.txt
    free -g >> scale_batch.node.txt

    start=\$(date +%s)
    blastn -task megablast -query ${queries} -db ${db_path} \\
        -evalue 1e-3 -max_target_seqs 5 -num_threads ${task.cpus} \\
        -outfmt '6 qseqid sgi sacc length pident mismatch gapopen qstart qend qlen sstart send slen sstrand evalue bitscore qcovhsp stitle staxids qseq sseq sseqid qcovs qframe sframe' \\
        -out scale_batch.bls
    end=\$(date +%s)

    echo "start=\${start}" > scale_batch.timing
    echo "end=\${end}" >> scale_batch.timing
    echo "seconds=\$(awk -v s=\${start} -v e=\${end} 'BEGIN { print e - s }')" >> scale_batch.timing
    """
}

process SCALE_SPLIT_CHUNK {
    tag "${chunk_id}"
    container params.blast_container
    publishDir params.outdir, mode: 'copy'
    cache false

    input:
    tuple val(chunk_id), path(chunk_fasta)
    val db_path
    val ready

    output:
    path "${chunk_id}.bls", emit: bls
    path "${chunk_id}.timing", emit: timing
    path "${chunk_id}.node.txt", emit: node

    script:
    """
    hostname > ${chunk_id}.node.txt
    nproc >> ${chunk_id}.node.txt
    free -g >> ${chunk_id}.node.txt

    start=\$(date +%s)
    blastn -task megablast -query ${chunk_fasta} -db ${db_path} \\
        -evalue 1e-3 -max_target_seqs 5 -num_threads ${task.cpus} \\
        -outfmt '6 qseqid sgi sacc length pident mismatch gapopen qstart qend qlen sstart send slen sstrand evalue bitscore qcovhsp stitle staxids qseq sseq sseqid qcovs qframe sframe' \\
        -out ${chunk_id}.bls
    end=\$(date +%s)

    echo "start=\${start}" > ${chunk_id}.timing
    echo "end=\${end}" >> ${chunk_id}.timing
    echo "seconds=\$(awk -v s=\${start} -v e=\${end} 'BEGIN { print e - s }')" >> ${chunk_id}.timing
    """
}
