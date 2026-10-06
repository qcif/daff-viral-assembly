#!/usr/bin/env nextflow
// Workflow (1): all queries in one FASTA, one blastn call on the whole node.

include { NODE_UP; BLAST_BATCH } from './modules.nf'

workflow {
    NODE_UP(params.blastn_db)
    BLAST_BATCH(file(params.queries), params.blastn_db, NODE_UP.out.done)

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
