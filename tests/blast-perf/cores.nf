#!/usr/bin/env nextflow
// Workflow (3): one node, one cold pass at 48 threads, then a sequential
// sweep of warm passes at increasing thread counts. All passes batch all
// 10 queries into one blastn call (the shape that won batch-vs-series.md);
// only -num_threads varies.

include {
    NODE_UP;
    BLAST_PASS as BLAST_PASS_COLD;
    BLAST_PASS as BLAST_PASS_WARM_1;
    BLAST_PASS as BLAST_PASS_WARM_2;
    BLAST_PASS as BLAST_PASS_WARM_3;
    BLAST_PASS as BLAST_PASS_WARM_4;
    BLAST_PASS as BLAST_PASS_WARM_5;
} from './modules.nf'

params.warm_threads = [2, 4, 8, 16, 48]

workflow {
    NODE_UP(params.blastn_db)
    BLAST_PASS_COLD(
        file(params.queries), params.blastn_db, 48, "cold", NODE_UP.out.done
    )

    // Five explicit sequential calls, not a loop: processes can't be
    // called from a `for`/closure in a workflow block. Each is gated on
    // the previous pass's `done` output, so they never overlap on this
    // one node. params.warm_threads must have exactly 5 entries.
    def (t1, t2, t3, t4, t5) = params.warm_threads
    BLAST_PASS_WARM_1(file(params.queries), params.blastn_db, t1, "warm", BLAST_PASS_COLD.out.done)
    BLAST_PASS_WARM_2(file(params.queries), params.blastn_db, t2, "warm", BLAST_PASS_WARM_1.out.done)
    BLAST_PASS_WARM_3(file(params.queries), params.blastn_db, t3, "warm", BLAST_PASS_WARM_2.out.done)
    BLAST_PASS_WARM_4(file(params.queries), params.blastn_db, t4, "warm", BLAST_PASS_WARM_3.out.done)
    BLAST_PASS_WARM_5(file(params.queries), params.blastn_db, t5, "warm", BLAST_PASS_WARM_4.out.done)

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
