# BLAST batch-vs-per-query performance test

Compares two ways of running MEGABLAST against `core_nt` on the `view`
Azure Batch pool: one task with all 10 queries in a single FASTA
(`batch.nf`), vs. one task per query run in parallel (`per_query.nf`).
See `task.md` for the full design rationale.

This is disposable: everything lives under this directory, the
production workflow is untouched, and the directory can be deleted once
you're done.

## Files

- `queries.fasta` / `queries.tsv` — 10 curated RefSeq viral genomes used
  as queries (committed, so both runs use byte-identical input).
  Regenerate with `./fetch_queries.sh` if needed.
- `modules.nf`, `batch.nf`, `per_query.nf` — the workflows.
- `blast-perf.config` — Azure config (real run, against `core_nt`).
- `smoke-test.config` — local Docker config (tiny test DB, for dry-running
  the wiring before spending anything on Azure).
- `summarise.py` — turns a run's outputs into `REPORT.md`.
- `run.sh` — drives all of the above.

## Local smoke test (no Azure cost)

```sh
nextflow run tests/blast-perf/batch.nf -c tests/blast-perf/smoke-test.config --outdir tests/blast-perf/results/smoke/batch
nextflow run tests/blast-perf/per_query.nf -c tests/blast-perf/smoke-test.config --outdir tests/blast-perf/results/smoke/per_query
```

## Real run (Azure, costs money)

Each workflow is its own invocation so the `view` pool can scale back to
0 nodes between them, rather than sharing one warm node. Run from the
repository root, and double check the `view` pool isn't already busy
before you start (`run.sh` will prompt you to confirm either way).

```sh
./tests/blast-perf/run.sh batch
# wait for the pool to scale back to 0 nodes (~15 min after its queue empties)
./tests/blast-perf/run.sh per_query <run_id>   # run_id printed by the batch step
./tests/blast-perf/run.sh report <run_id>
```

The report lands at `results/<run_id>/REPORT.md`.
