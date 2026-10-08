# BLAST performance tests

Three disposable Azure Batch performance tests for MEGABLAST against
`core_nt` on the `view` pool. Everything lives under this directory, the
production workflow is untouched, and the directory can be deleted once
you're done.

## 1. Batch vs. per-query

Compares one task with all 10 queries in a single FASTA (`batch.nf`) vs.
one task per query run in parallel (`per_query.nf`). See
`batch-vs-series.md` for the full design rationale.

### Files

- `queries.fasta` / `queries.tsv` — 10 curated RefSeq viral genomes used
  as queries (committed, so both runs use byte-identical input).
  Regenerate with `./fetch_queries.sh` if needed.
- `modules.nf`, `batch.nf`, `per_query.nf` — the workflows.
- `blast-perf.config` — Azure config (real run, against `core_nt`).
- `smoke-test.config` — local Docker config (tiny test DB, for dry-running
  the wiring before spending anything on Azure).
- `summarise.py` — turns a run's outputs into `REPORT.md`.

### Local smoke test (no Azure cost)

```sh
nextflow run tests/blast-perf/batch.nf -c tests/blast-perf/smoke-test.config --outdir tests/blast-perf/results/smoke/batch
nextflow run tests/blast-perf/per_query.nf -c tests/blast-perf/smoke-test.config --outdir tests/blast-perf/results/smoke/per_query
```

### Real run (Azure, costs money)

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

## 2. Thread-count sweep

Measures MEGABLAST speed (batched shape) across `-num_threads` 2, 4, 8,
16 and 48, once with `core_nt` cold and once per thread count with it
warm (cached). See `cores-allocated.md` for the full design rationale.

### Files

- `modules.nf`'s `BLAST_PASS` process, `cores.nf` — the workflow (one
  cold pass then a sequential warm sweep, all on one node).
- `blast-perf.config` / `smoke-test.config` — same config files as above;
  each has a `BLAST_PASS_.*` block.
- `summarise_cores.py` — turns a run's outputs into `REPORT.md`.

### Local smoke test (no Azure cost)

```sh
nextflow run tests/blast-perf/cores.nf -c tests/blast-perf/smoke-test.config --outdir tests/blast-perf/results/smoke/cores
```

### Real run (Azure, costs money)

One invocation covers the whole sweep (cold pass + all five warm passes),
so there's no scale-down wait within this test.

```sh
./tests/blast-perf/run.sh cores
./tests/blast-perf/run.sh cores-report <run_id>   # run_id printed above
```

The report lands at `results/<run_id>/REPORT.md`.

## 3. Scale test (~30k queries)

The two tests above both use only 10 queries, so they mostly measure the
cost of scanning `core_nt` once, not the cost of aligning many queries.
This test repeats the batch-vs-split comparison at production's real
query volume: all ~30k queries in one task (`scale_batch.nf`, 48 cores)
vs. the query file split into 10 chunks run in parallel (`scale_split.nf`,
16 cores each, 3 at a time). See `scale-test.md` for the full design
rationale.

### Files

- `fetch_queries_scale.sh` — fetches ~30,000 viral nucleotide sequences
  (9-11kb) from NCBI into `scale_queries/queries_30k.fasta`. Not
  committed (hundreds of MB; see `.gitignore`) — regenerate as needed.
  Takes several minutes (NCBI's unauthenticated rate limit).
- `modules.nf`'s `SCALE_BATCH`/`SCALE_SPLIT_CHUNK` processes,
  `scale_batch.nf`, `scale_split.nf` — the workflows.
- `blast-perf.config` / `smoke-test.config` — same config files as
  above; each has `SCALE_BATCH`/`SCALE_SPLIT_CHUNK` blocks.
- `summarise_scale.py` — turns a run's outputs into `REPORT.md`.

### Local smoke test (no Azure cost, no fetch needed)

Uses the small committed `queries.fasta` (10 queries → 10 single-query
chunks), just to prove the wiring:

```sh
nextflow run tests/blast-perf/scale_batch.nf -c tests/blast-perf/smoke-test.config --outdir tests/blast-perf/results/smoke/scale_batch
nextflow run tests/blast-perf/scale_split.nf -c tests/blast-perf/smoke-test.config --outdir tests/blast-perf/results/smoke/scale_split
```

### Real run (Azure, costs money — more than the other two tests)

```sh
./tests/blast-perf/fetch_queries_scale.sh   # if scale_queries/queries_30k.fasta doesn't exist yet
./tests/blast-perf/run.sh scale_batch
# wait for the pool to scale back to 0 nodes (~15 min after its queue empties)
./tests/blast-perf/run.sh scale_split <run_id>   # run_id printed by the scale_batch step
./tests/blast-perf/run.sh scale-report <run_id>
```

The report lands at `results/<run_id>/REPORT.md`. Query volume is ~3,000x
the other two tests, so time and cost are far less predictable — see the
estimate and caveats in `scale-test.md`.
