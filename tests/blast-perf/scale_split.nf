#!/usr/bin/env nextflow
// Scale test, shape (b): the ~30k queries split into params.n_chunks (10)
// roughly-equal batches, each run as its own blastn call with 16 cores,
// up to params.chunk_max_forks (3) running at once — mirrors the 16C/3x
// slot-packing shape from batch-vs-series.md, but as a small fixed number
// of batches rather than one task per query. Companion to scale_batch.nf
// — run each on its own fresh node (see README.md).

include { NODE_UP; SCALE_SPLIT_CHUNK } from './modules.nf'

params.n_chunks = 10

workflow {
    NODE_UP(params.blastn_db)

    // Chunk size computed from the actual query count so n_chunks comes
    // out exact regardless of how many records fetch_queries_scale.sh
    // ended up keeping. `grep -c` runs locally (where `nextflow run` is
    // launched), not on Azure — same as any other up-front Groovy in a
    // workflow block.
    def n_queries = "grep -c ^> ${params.queries}".execute().text.trim() as int
    def n_chunks = params.n_chunks as int
    def chunk_size = Math.ceil(n_queries / n_chunks.doubleValue()) as int
    log.info "scale_split: ${n_queries} queries / ${params.n_chunks} chunks" +
        " -> ${chunk_size} queries per chunk"

    ch_chunks = Channel.fromPath(params.queries)
        .splitFasta(by: chunk_size, file: true)
        .map { f -> tuple(f.baseName, f) }
    SCALE_SPLIT_CHUNK(ch_chunks, params.blastn_db, NODE_UP.out.done)

    workflow.onComplete = {
        def timings = [
            start: workflow.start.toString(),
            complete: workflow.complete.toString(),
            duration: workflow.duration.toString(),
            success: workflow.success,
        ]
        file("${params.outdir}/timings.json")
            .text = groovy.json.JsonOutput.toJson(timings)
    }
}
