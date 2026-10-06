#!/usr/bin/env nextflow
// Workflow (2): one task per query, run in parallel.

include { NODE_UP; BLAST_SINGLE } from './modules.nf'

workflow {
    NODE_UP(params.blastn_db)
    ch_queries = Channel.fromPath(params.queries)
        .splitFasta(by: 1, file: true)
        .map { f -> tuple(f.baseName, f) }
    BLAST_SINGLE(ch_queries, params.blastn_db, NODE_UP.out.done)

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
