# Plan: MEGABLAST vs core_nt performance test on Azure Batch

## Goal

Find out which is faster on the `view` pool: one batched BLAST task, or one
task per query. Both run MEGABLAST against `core_nt` with the same 10 viral
queries (about 10kb each), and each is tuned for the node rather than copied
from production settings.

| Workflow | Shape | Resources |
|---|---|---|
| `batch.nf` | 1 task: all 10 queries in one FASTA → one `blastn` call | 48 CPUs (whole node) |
| `per_query.nf` | 10 tasks: one query each, run in parallel | 16 CPUs each (Batch fits 3 at a time) |

Each workflow runs **once, with core_nt evicted from the page cache** before
BLAST starts. Both runs then read the database from NVMe, as a production run
on a fresh node does. They share one node, back-to-back, so there's no need to
wait for scale-down or re-stage 830GB of refdata between runs.

The deliverable is `tests/blast-perf/REPORT.md`, which reports the walltime of
each workflow and states which one is faster.

This is disposable. Everything goes under `tests/blast-perf/`, the production
workflow is not changed, and the directory can be deleted when the test is
done.

## Infrastructure facts this relies on

- Pool `view`: one `Standard_L48as_v3` node (48 vCPU, 384GB RAM, NVMe),
  `taskSlotsPerNode=48`. It autoscales 0→1 when tasks are pending, and scales
  back to 0 about 15 minutes after the queue empties
  (`deploy/azure/batch-helpers.sh`).
- The start task stages about 830GB of refdata to `/mnt/nvme/refdata/`,
  including `core_nt` at `/mnt/nvme/refdata/core_nt/core_nt` (269GB).
- Container: `quay.io/biocontainers/blast:2.16.0--h66d330f_4`, the same image
  production uses.

## Directory layout

```
tests/blast-perf/
├── queries.fasta          # 10 curated viral genomes (committed)
├── queries.tsv            # accession, name, family, length
├── fetch_queries.sh       # regenerates queries.fasta from NCBI
├── modules.nf             # EVICT_DB, BLAST_BATCH, BLAST_SINGLE
├── node_up.nf             # provisions the node (one EVICT_DB task)
├── batch.nf               # workflow (1)
├── per_query.nf           # workflow (2)
├── blast-perf.config      # Azure config for all three workflows
├── run.sh                 # node_up → batch → per_query, collects results
├── summarise.py           # traces → REPORT.md
├── results/<run_id>/      # trace.txt, timings.json, logs per run
└── REPORT.md              # deliverable
```

## 1. Curate the queries

Use RefSeq complete genomes, 9–11kb each. Pick a mix of plant viruses (the
VIEW use case) and well-sampled animal viruses. Some queries will then return
many core_nt hits and others fewer, which is closer to real contigs than 10
near-identical queries would be.

| # | Accession | Virus | Family | ~Length |
|---|---|---|---|---|
| 1 | NC_001616 | Potato virus Y | Potyviridae | 9.7kb |
| 2 | NC_001445 | Plum pox virus | Potyviridae | 9.7kb |
| 3 | NC_002509 | Turnip mosaic virus | Potyviridae | 9.8kb |
| 4 | NC_003224 | Zucchini yellow mosaic virus | Potyviridae | 9.6kb |
| 5 | NC_002634 | Soybean mosaic virus | Potyviridae | 9.6kb |
| 6 | NC_001886 | Wheat streak mosaic virus | Potyviridae | 9.4kb |
| 7 | NC_001477 | Dengue virus 1 | Flaviviridae | 10.7kb |
| 8 | NC_012532 | Zika virus | Flaviviridae | 10.8kb |
| 9 | NC_009942 | West Nile virus | Flaviviridae | 11.0kb |
| 10 | NC_002031 | Yellow fever virus | Flaviviridae | 10.9kb |

`fetch_queries.sh` does the following:

1. Fetches each accession with `efetch -db nucleotide -format fasta` (or with
   E-utilities over `curl` if entrez-direct is not installed).
2. Rewrites each header to `>q01_NC_001616`, `>q02_...` and so on, so the IDs
   are short, unique and sort in order.
3. Checks that each length is between 9,000 and 11,500bp and fails otherwise.
   Replace any accession that is out of range or has become obsolete.
4. Writes `queries.tsv`.

Commit `queries.fasta` so both runs use byte-identical input.

## 2. Evicting core_nt from memory

### Why this is needed even on a fresh node

The start task writes core_nt to NVMe with azcopy. Those writes pass through
the page cache, so a "cold" node can still hold a large, unpredictable share
of core_nt in RAM by the time the first task runs. A fresh node doesn't give a
controlled cold cache. Explicit eviction does, and it also lets both runs
share one node.

### Method: `posix_fadvise(POSIX_FADV_DONTNEED)` per DB file

`EVICT_DB` walks `/mnt/nvme/refdata/core_nt/` and calls
`os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_DONTNEED)` on every file. This asks
the kernel to drop those files' clean pages from the page cache.

- **No privileges needed.** It works on a read-only fd, so the existing
  `:ro` bind mount is enough: no root, no `--privileged` container, and no
  changes to the pool. The page cache belongs to the host kernel, so eviction
  from inside the container affects the whole node.
- **Targeted.** Only core_nt is dropped. It doesn't disturb anything else on
  the node.
- **Verified.** The task runs `fincore` (util-linux) on the DB files before
  and after, and records `free -g`. Its log (`evict.log`) gives the cached
  bytes before and after. If more than 1% of core_nt is still resident after
  eviction, the task fails.
- **Container.** Use a fully qualified image with Python 3 and util-linux,
  e.g. `docker.io/library/python:3.12-slim` (Debian base, includes `fincore`).
  Batch rejects unqualified names (see `azure-batch-handover.md`).

**Fallback** if `fadvise` leaves pages resident (it skips pages that are dirty
or mapped by a running process, but nothing else should be using core_nt):
`sync; echo 1 > /proc/sys/vm/drop_caches`. This drops the whole node's page
cache. It needs `containerOptions '--privileged'` on `EVICT_DB` alone. Use it
only if the `fincore` check fails.

### Ordering

Each workflow runs `EVICT_DB` first and gates BLAST on its output, so no BLAST
task starts until eviction is done and verified:

```groovy
EVICT_DB(params.blastn_db)
BLAST_BATCH(file(params.queries), params.blastn_db, EVICT_DB.out.done)
```

`BLAST_SINGLE` takes `EVICT_DB.out.done` as a value channel too, so all 10
tasks wait on the same eviction.

## 3. Shared BLAST command

Both workflows use the search settings from `modules/blast/blastn/main.nf`, so
the hit workload is representative of production:

```
blastn -task megablast -query <q> -db ${blastn_db} \
    -evalue 1e-3 -max_target_seqs 5 -num_threads ${task.cpus} \
    -outfmt '6 qseqid sgi sacc length pident ... qframe sframe'
```

It differs from the production module in three ways:

- **The DB is a plain string (`val`).** It's passed straight to `-db`,
  matching production's `database_mode = "mounted"` behaviour under the
  `azure` profile.
- **`-num_threads ${task.cpus}`** instead of `${params.blast_threads}`, so
  BLAST uses every core the task was given.
- **In-task timing.** The process script brackets `blastn` with
  `date +%s.%N` and writes `<id>.timing` (start, end and seconds). It also
  writes `node.txt` with `hostname`, `nproc` and `free -g`.

All processes set `cache false` so that a re-run never resumes.

`-mt_mode` stays at the default `0` (threads split the database). Mode `1`
threads by query, and with only 10 queries it would cap the 48-CPU batch run
at 10 busy threads.

## 4. Workflows

The processes live in `modules.nf`. Each workflow `include`s what it needs.

### `batch.nf`

```groovy
workflow {
    EVICT_DB(params.blastn_db)
    BLAST_BATCH(file(params.queries), params.blastn_db, EVICT_DB.out.done)
}
```

### `per_query.nf`

```groovy
workflow {
    EVICT_DB(params.blastn_db)
    ch_queries = Channel.fromPath(params.queries)
        .splitFasta(by: 1, file: true)
        .map { f -> tuple(f.baseName, f) }
    BLAST_SINGLE(ch_queries, params.blastn_db, EVICT_DB.out.done)
}
```

### `node_up.nf`

This is a single `EVICT_DB` task. Its only job is to make autoscale provision
the node and run the start task, so that neither measured workflow pays for
node start-up. Its trace is kept, and the report shows provisioning time
separately.

Both measured workflows write `timings.json` from `workflow.onComplete`
(`workflow.start`, `workflow.complete`, `workflow.duration`, `workflow.success`)
into `${params.outdir}`. Every run is launched with
`-with-trace -with-report`, using these trace fields:
`task_id,name,tag,status,exit,submit,start,complete,duration,realtime,%cpu,peak_rss,rchar`.

## 5. Config: `blast-perf.config`

Reuse `conf/azure.config` for the `azure {}` block, the executor, the queue
(`view`) and `containerOptions` via `includeConfig`, then set only what the
test needs:

```groovy
includeConfig '../../conf/azure.config'

workDir = 'az://workdata/blast-perf'

params {
    queries = "${projectDir}/queries.fasta"
    blastn_db = '/mnt/nvme/refdata/core_nt/core_nt'
}

process {
    withName: 'EVICT_DB' {
        container = 'docker.io/library/python:3.12-slim'
        cpus = 1
        memory = 2.GB
        time = 30.m
    }
    withName: 'BLAST_BATCH' {
        cpus = 48
        memory = 300.GB
        time = 4.h
    }
    withName: 'BLAST_SINGLE' {
        cpus = 16
        memory = 100.GB
        time = 4.h
    }
}
```

### Parallelism without `maxForks`

Nextflow doesn't cap concurrency here. It submits all 10 `BLAST_SINGLE` tasks
to Batch at once. Azure Batch then packs them onto the node by **task slots**:
Nextflow gives each task a slot count based on its share of the VM's CPUs
**and memory**, whichever is larger. The pool has `taskSlotsPerNode=48`.

- `BLAST_SINGLE`: 16/48 of the CPUs and 100/384 of the memory → 16 slots →
  **3 run at once**, the rest queue in Batch → 4 waves (3 + 3 + 3 + 1).
- `BLAST_BATCH`: 48 CPUs → all 48 slots.

So yes, it should come out as 3 × 16C. The catch is that **memory also drives
the slot count**. A `BLAST_SINGLE` memory request above 128GB (a third of
384GB) would raise each task to more than 16 slots, leaving room for only 2 at
once. Keep it at or below 128GB. `summarise.py` confirms this from the trace:
it counts the maximum number of `BLAST_SINGLE` tasks with overlapping
start–complete intervals and flags anything other than 3.

The new process names don't collide with any `withName` key in
`azure.config`, so the replace-not-merge problem described in the handover
can't drop the container. The BLAST container is also declared in each
process as a second safeguard.

Launch from the repo root:

```sh
nextflow run tests/blast-perf/batch.nf -c tests/blast-perf/blast-perf.config ...
```

## 6. Run protocol

`run.sh` runs everything back-to-back on one node:

1. **Pre-flight.** Source `.env.azure` and `deploy/azure/batch-helpers.sh`.
   Check that `view` has no active jobs (`az_jobs_list`), because this is the
   production pool. Abort if it's busy.
2. **`node_up.nf`.** This provisions the node (VM boot + start task + refdata
   staging). Its trace is recorded for the report but isn't part of the
   comparison.
3. **`batch.nf`.** Evict, then BLAST. Launch it immediately after step 2, so
   the node is still inside its 15-minute keep-warm window.
4. **`per_query.nf`.** Evict, then BLAST. Launch it immediately after step 3.
5. **Collect** each run into `results/<workflow>_<timestamp>/` (trace, report,
   timings.json, `evict.log`, `.bls`, `.timing`, `node.txt`).
6. **Tear-down.** Nothing to do: autoscale returns the pool to 0 after about
   15 minutes idle. Delete the `blast-perf/` prefix in the `workdata`
   container, or leave it for the 14-day lifecycle policy.

**Same-node check:** `summarise.py` confirms that every task across both
measured runs reports the same hostname. If the node was recycled in between,
the second run paid for provisioning, and the report flags it.

**Remaining asymmetry:** `batch.nf` pulls the BLAST container image and
`per_query.nf` reuses it. That costs seconds against a run of many minutes,
and it falls in the `submit → start` gap, which is reported separately. To
remove it, have `node_up.nf` also run a trivial task in the BLAST image.

## 7. Correctness check

Speed only counts if both shapes give the same answer. `summarise.py` sorts the
batch `.bls` file and the concatenated per-query `.bls` files and diffs them.
It reports identical results, or the number of differing rows per query.
Megablast batching can occasionally change hits at the `max_target_seqs`
cut-off, so a small difference is worth recording rather than treating as a
failure.

## 8. Report: `summarise.py` → `REPORT.md`

`summarise.py` uses only the standard library (`csv`, `json`, `datetime`),
runs from `./venv`, and is checked with flake8.

| Metric | Source | Meaning |
|---|---|---|
| **Workflow walltime** | `timings.json` duration | Headline number: launch → done, node already up |
| **BLAST span** | first `blastn` start → last `blastn` end (`.timing`) | Pure compute comparison, excluding eviction |
| Eviction | `evict.log` | Time taken, and core_nt resident before → after |
| Σ blastn seconds × cpus | `.timing` + config | CPU-seconds of BLAST work |
| Per-task Azure overhead | (submit→start) + (start→blastn start) | Scheduling + container + staging cost of splitting |
| Max concurrent `BLAST_SINGLE` | trace intervals | Confirms the expected 3 × 16C |
| Per-query blastn time + wave | `.timing` (per_query only) | Shows whether one query dominates the last wave |
| Peak RSS, %CPU | trace | Whether threads were actually used |
| Est. node cost | walltime × $7.50/h | |

The report contains:

- A summary table with one row per workflow.
- A per-query table for `per_query`.
- The eviction, concurrency, same-node and correctness check results.
- Node provisioning time (from `node_up.nf`) and the hostname, for context.
- A short **verdict** naming the faster shape on this infrastructure and by
  what margin on both walltime and BLAST span.

The verdict also has to say plainly that each figure comes from **a single
run**: there's no measure of run-to-run spread, so a small margin (under about
10–15%) should be treated as inconclusive.

## Implementation order

1. `fetch_queries.sh` → `queries.fasta` and `queries.tsv` (verify lengths).
2. `modules.nf`, the three workflows and `blast-perf.config`. Smoke-test
   locally with a docker override against `tests/blast_test_db/Wsmv.fa`, with
   CPUs scaled down. Query 6 (WSMV) should hit. `EVICT_DB` should show the
   test DB resident before eviction and gone after, which proves the
   `fadvise` + `fincore` approach before any Azure spend.
3. `summarise.py`, run against the local smoke-test outputs; flake8 clean.
4. `run.sh`.
5. Azure: `run.sh` (node_up → batch → per_query).
6. Generate `REPORT.md` and review it.

## Cost and time estimate

Node start + refdata staging is likely 15–40 minutes. Each measured run is
eviction (seconds) + BLAST from NVMe (likely 5–20 minutes). Add about
15 minutes of keep-warm at the end. That is roughly 1–1.5 node-hours in
total, about $8–12 at $7.50/h, with no idle wait between runs.

## Follow-ups (outside this test)

1. **`-num_threads ${params.blast_threads}` (default 2) against `cpus = 16`.**
   This is in `modules/blast/blastn/main.nf`. Switch to `${task.cpus}`.
2. **Resources.** Retune the production MEGABLAST `cpus` (and its `maxForks`,
   which may then be unnecessary given slot packing) based on the result.
