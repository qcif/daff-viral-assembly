#!/usr/bin/env nextflow
// Scale test, shape (a): one blastn call over all ~30k queries at once,
// using the whole node (48 cores). Companion to scale_split.nf — run each
// on its own fresh node (see README.md), not back-to-back on one node,
// so neither pays for the other's warm cache.

include { NODE_UP; SCALE_BATCH } from './modules.nf'

workflow {
    NODE_UP(params.blastn_db)
    SCALE_BATCH(file(params.queries), params.blastn_db, NODE_UP.out.done)

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
