# Plan: MEGABLAST at production query volume (~30k) on Azure Batch

## Goal

The batch-vs-per-query test (`batch-vs-series.md`) and the thread-count
sweep (`cores-allocated.md`) both used only 10 queries, so almost all of
the measured time was the cost of scanning `core_nt` once, not the cost
of aligning many queries. Neither result can be assumed to hold at
production's real query volume (~30,000 contigs per run). This test
re-runs the batch-vs-split comparison at that scale, to see whether the
earlier verdict ("batch wins") still holds once per-query cost actually
dominates.

| Run | Shape | Resources |
|---|---|---|
| `scale_batch.nf` | 1 task: all ~30k queries in one FASTA → one `blastn` call | 48 CPUs (whole node) |
| `scale_split.nf` | 10 tasks: query file split into 10 equal chunks, run in parallel | 16 CPUs each, 3 at a time (`maxForks = 3`) |

Each workflow runs once, on its own fresh node (pool scaled back to 0
between them), exactly as `batch.nf`/`per_query.nf` did.

This reuses the existing harness in `tests/blast-perf/` and stays
disposable.

## Why split into 10 fixed chunks, not one task per query

`per_query.nf` (the first test) ran one task per query — fine at 10
queries, but 30,000 Azure Batch tasks is a different scheduling regime
(submission overhead, Batch API limits, pool queue depth) that isn't what
this test is trying to measure. Splitting into a small, fixed number of
batches is closer to how this would actually be deployed if per-query
parallelism won, and keeps the task count comparable to the node's real
concurrency (3 tasks fit at 16 CPUs each on a 48-CPU node, as established
in `batch-vs-series.md`).

## Query set

~30,000 viral nucleotide sequences, 9,000-11,000bp, fetched from NCBI
nuccore (not just RefSeq — RefSeq alone only has ~900 records in this
length window) via `fetch_queries_scale.sh`:

```sh
./tests/blast-perf/fetch_queries_scale.sh
```

This is **not curated** the way the 10-query set was (no attempt to
balance virus families or exclude near-duplicate strains) — at this
volume the goal is realistic total data volume and per-query cost, not a
controlled set of expected hits. It takes several minutes (NCBI's
unauthenticated rate limit is 3 requests/sec, and there's no API key
configured for this project) and writes:

- `scale_queries/queries_30k.fasta` — the sequences, ~250-350MB.
- `scale_queries/queries_30k.n` — exact record count (filtering during
  fetch can land slightly under or over 30,000).

Neither is committed (`.gitignore`): regenerate as needed. The exact
records returned can drift between runs as NCBI's database changes,
which is fine for a throughput test but means this test doesn't have the
10-query test's exact-correctness guarantee across runs — only within a
single pair of runs (batch vs. split both read the one fetched file).

## Workflows

`scale_split.nf` computes the chunk size itself (`ceil(n_queries / 10)`
via `grep -c '^>'` on the query file, run locally before launching, same
as any other up-front Groovy in a workflow block) so splitting into
exactly 10 chunks doesn't depend on the fetch landing on exactly 30,000
records.

```groovy
// scale_batch.nf
NODE_UP(params.blastn_db)
SCALE_BATCH(file(params.queries), params.blastn_db, NODE_UP.out.done)

// scale_split.nf
NODE_UP(params.blastn_db)
def chunk_size = Math.ceil(n_queries / 10.0) as int
ch_chunks = Channel.fromPath(params.queries)
    .splitFasta(by: chunk_size, file: true)
    .map { f -> tuple(f.baseName, f) }
SCALE_SPLIT_CHUNK(ch_chunks, params.blastn_db, NODE_UP.out.done)
```

Both processes (`SCALE_BATCH`, `SCALE_SPLIT_CHUNK` in `modules.nf`) are
the same shape as `BLAST_BATCH`/`BLAST_SINGLE` from the first test: same
`blastn` command and `-outfmt`, in-task timing via `date +%s`, `hostname`/
`nproc`/`free -g` into `node.txt`. `cpus`/`memory`/`maxForks` come from
`blast-perf.config`:

| Process | cpus | memory | maxForks |
|---|---|---|---|
| `SCALE_BATCH` | 48 | 300.GB | 1 |
| `SCALE_SPLIT_CHUNK` | 16 | 100.GB | 3 |

`time` is raised to `6.h` for both (vs. `4.h` in the first test) since
30,000 queries is roughly 3,000x the first test's volume.

## Run protocol (`run.sh`)

```sh
./tests/blast-perf/run.sh scale_batch
# wait for the pool to scale back to 0 nodes (~15 min after its queue empties)
./tests/blast-perf/run.sh scale_split <run_id>
./tests/blast-perf/run.sh scale-report <run_id>
```

`run.sh` passes `--queries tests/blast-perf/scale_queries/queries_30k.fasta`
for both modes (overriding `blast-perf.config`'s default, which still
points `batch.nf`/`per_query.nf`/`cores.nf` at the small 10-query set).

## Report: `summarise_scale.py` → `REPORT.md`

Same shape as `summarise.py`'s report (imports its helpers), adapted for
one batch task vs. 10 split chunks instead of 10 per-query tasks:
workflow walltime, BLAST span, summed CPU-seconds, per-chunk timing,
max-concurrent-chunks check (expect 3), correctness (sorted `.bls` diff
between the batch output and the concatenated chunk outputs), node
provisioning per run, verdict, estimated cost.

## Implementation order

1. `fetch_queries_scale.sh` → `scale_queries/queries_30k.fasta` (run in
   background; takes several minutes).
2. `SCALE_BATCH`/`SCALE_SPLIT_CHUNK` in `modules.nf`, `scale_batch.nf`,
   `scale_split.nf`, config blocks in `blast-perf.config` and
   `smoke-test.config`.
3. Smoke-test locally against the existing 10-query `queries.fasta`
   (`scale_split.nf --queries queries.fasta` then gives 10 single-query
   chunks) — proves the chunk-size computation and wiring without
   waiting on the real fetch.
4. `summarise_scale.py` against the smoke outputs; flake8.
5. `run.sh` modes; README section.
6. Once the fetch completes: Azure `run.sh scale_batch`, wait for
   scale-down, `run.sh scale_split`, `run.sh scale-report`.

## Cost and time estimate

Per invocation: ~15-20 minutes node start-up + refdata staging, plus
~15 minutes keep-warm at the end. The BLAST work itself is the open
question this test exists to answer, but as a rough floor: the first
test's 10 queries took ~2-5 minutes depending on shape, so 30,000 queries
(3,000x the volume, though not purely linear — see `cores-allocated.md`
on diminishing per-thread returns) could run from tens of minutes to
several hours. Estimated **$10-20 per workflow**, but this is a wide
range until the first run narrows it down — if a run is clearly headed
past a few hours, consider killing it rather than waiting it out.

## Follow-ups (outside this test)

1. If `scale_split.nf` wins decisively, it's worth also testing
   `-mt_mode 1` (threads split by query, not by database) — the NCBI
   recommendation for large query counts — since this test still uses the
   default `-mt_mode 0` throughout, same as the first two tests.
2. These are still RefSeq-style complete genomes, not real contigs. If
   either shape's result is close, re-test with actual production contig
   output (varied lengths, more hits per query, host sequence mixed in)
   before changing production settings.
