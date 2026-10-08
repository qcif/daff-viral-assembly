# Plan: MEGABLAST thread count on Azure Batch

## Goal

Find out how much MEGABLAST against `core_nt` on the `view` pool actually
benefits from more threads, cold and warm. Production currently requests
`blast_threads: 2` (`params/default_params.yml`) but `cpus` up to 48
depending on the process. If few threads are nearly as fast once the
database is cached, production can pack many more BLAST tasks onto a node.

| Pass | Cache state | Threads tried |
|---|---|---|
| Cold | first read of `core_nt` on a fresh node | 48 only |
| Warm | `core_nt` fully cached (384GB node, DB is 269GB) | 2, 4, 8, 16, 48 |

All passes use the batched shape (all 10 queries in one `blastn` call),
which won the batch-vs-per-query test (`batch-vs-series.md`,
`results/20261006_081816/REPORT.md`). Only `-num_threads` changes between
passes.

This reuses the existing harness in `tests/blast-perf/` and stays
disposable, like the first test.

## What the last run tells us

The raw per-task trace from `20261006_081816` (pulled from
`az://workdata/blast-perf/`, since `.command.trace` isn't published by
default) gives CPU time and bytes actually read from disk, not just
Nextflow's `%cpu` average:

| Task | Cache | Threads | Time | Avg %CPU | Disk read |
|---|---|---|---|---|---|
| `BLAST_BATCH`, 10 queries | cold | 48 | 135s | 453% | 264GB |
| `BLAST_SINGLE` wave 1 (q5, q9, q10) | cold | 16 each | ~180s | ~165% | ~88GB each |
| `BLAST_SINGLE` later waves | warm | 16 each | ~26s | ~1,100% | <1GB |

**`%cpu` in Nextflow's trace is a whole-task average (CPU-seconds used ÷
elapsed time), not a peak** — Nextflow has no peak-CPU sampling, only
`peak_rss`/`peak_vmem` for memory. A single number can't distinguish "used
all cores the whole time" from "burst to 48 cores then idled waiting on
disk", which matters a lot for the cold pass.

Reading those two prior numbers differently than the first draft of this
plan did:

- **Cold is disk-bound, not CPU-bound.** 264GB in 135s is about 2GB/s from
  NVMe. At 453% average CPU (~4.5 of 48 cores), most of the 48 threads
  were idle most of the time. More threads only help the cold pass if
  `-mt_mode 0` (the default) actually issues more concurrent reads — that's
  the open question, not raw compute.
- **Batching is cheap once warm.** The 48-thread batch call used about
  **611 CPU-seconds total** (453% × 135s) for all 10 queries combined.
  A single warm 16-thread per-query task used about 290 CPU-seconds
  — *per query* — because each query re-scans the whole database. So the
  warm, batched cost for 10 queries is roughly 2x one per-query task, not
  10x. At 4 threads that's roughly 611 / 4 ≈ 150s, not the ~12 minutes the
  first draft of this plan estimated.
- **The two workloads may be close at low thread counts.** With ~150s of
  compute at 4 threads and a ~135s disk floor on the cold pass, 4 and 48
  cores could come out similar, especially cold. A single 4-vs-48
  comparison risks a non-finding; a sweep shows the actual curve.

## Design

### One workflow, one invocation, a sweep of warm passes

New workflow `cores.nf`, run **once**, on one fresh node:

```groovy
workflow {
    NODE_UP(params.blastn_db)
    BLAST_PASS_COLD(file(params.queries), params.blastn_db, 48, "cold", NODE_UP.out.done)

    ch_threads = Channel.of(2, 4, 8, 16, 48)
    BLAST_PASS_WARM(file(params.queries), params.blastn_db, ch_threads, "warm", BLAST_PASS_COLD.out.done)
}
```

- **Cold pass:** one call at 48 threads, the first read of `core_nt` on
  this node. Only one cold data point is possible per node, so it isn't
  swept.
- **Warm sweep:** five calls at increasing thread counts, run
  **sequentially** (each gated on the previous one's `done` output) so
  they never share the node's CPUs. The database is cached after the cold
  pass, so all five are warm.

`NODE_UP` stays as the first task. It doesn't touch `core_nt`'s data files
(`blastdbcmd -info` only reads the index header), so the cold pass is still
cold. **Check this when building:** if `NODE_UP` turns out to warm the
cache noticeably, drop the `blastdbcmd` line from it for this test.

### Process changes (`modules.nf`)

Add one process, `BLAST_PASS`, aliased twice on include
(`include { BLAST_PASS as BLAST_PASS_COLD; BLAST_PASS as BLAST_PASS_WARM }`).
It's the same shape as `BLAST_BATCH` except:

- Takes `val threads` and `val pass` (`"cold"` / `"warm"`) inputs, and
  uses `-num_threads ${threads}` instead of `${task.cpus}` — the task's
  `cpus` stays fixed at 48 throughout (see below), only the `blastn` flag
  changes. Names outputs `${pass}_${threads}.bls` /
  `${pass}_${threads}.timing` / `${pass}_${threads}.node.txt`.
- Emits `val true, emit: done` so passes can be chained sequentially.
- **CPU/disk sampling.** Runs a background loop for the duration of
  `blastn`, once per second, appending to `${pass}_${threads}.sample`:
  `blastn`'s CPU time from `/proc/<pid>/stat` (fields 14+15, utime+stime)
  and total disk sectors read from `/proc/diskstats`, both as running
  counters (deltas computed in `summarise_cores.py`, not in the task).
  This is what actually answers "did more threads get used, or did it
  just wait on disk" — the single end-of-task average can't.
  ```
  blastn ... &
  bpid=$!
  while kill -0 \$bpid 2>/dev/null; do
      t=\$(date +%s)
      cpu=\$(awk '{print \$14+\$15}' /proc/\$bpid/stat 2>/dev/null || echo 0)
      sectors=\$(awk '\$3=="nvme0n1"{print \$6}' /proc/diskstats)
      echo "\$t \$cpu \$sectors" >> ${pass}_${threads}.sample
      sleep 1
  done
  wait \$bpid
  ```
  Adjust the disk device name to whatever `lsblk` shows for the refdata
  mount on pool `view` — confirm this during smoke-testing (the local
  smoke test has no NVMe device, so guard with `|| echo 0`, which also
  means the sampling loop only gives a real disk signal on Azure).
- `node.txt` also records `cat /sys/fs/cgroup/cpu.max` and
  `cat /sys/fs/cgroup/memory.max`. In the last run, `nproc` reported 48
  inside a 16-CPU task, so it isn't known whether Batch enforces the CPU
  request at all, or only `-num_threads` controls what BLAST actually
  uses. (cgroup v2 only; record `n/a` if the files are missing.)

Leave `BLAST_BATCH` and `BLAST_SINGLE` alone, so the first test can still
be re-run.

### Holding everything else constant

| Setting | Value | Why |
|---|---|---|
| `task.cpus` | 48, fixed for every pass | Keeps Batch's slot/scheduling behaviour identical across passes — only `-num_threads` (a `blastn` flag, not a Nextflow resource) varies |
| `memory` | 300.GB | Matches the first test; large enough that nothing gets evicted under memory pressure between passes |
| `maxForks` | 1 | One BLAST task on the node at a time; passes still run sequentially via the `done` chain regardless |
| `-mt_mode` | default (0) | Same as production and as the first test |
| Queries | `queries.fasta` | Same 10 queries as the first test |

### Config (`blast-perf.config`)

Add a `withName: 'BLAST_PASS_.*'` block: `cpus = 48`, `memory = 300.GB`,
`time = 1.h`, `maxForks = 1`. Restate the container, as the other blocks
do (see the replace-not-merge note in `blast-perf.config`).

`smoke-test.config` needs the same block scaled down (`cpus = 2`,
`memory = 1.GB`), and the thread sweep in `cores.nf` should shrink to
`Channel.of(1, 2)` under a `params.smoke` flag, or just accept that the
smoke test only proves the wiring, not the thread scaling (the tiny test
DB is too small for thread count to matter anyway).

## Run protocol (`run.sh`)

Add one mode:

```sh
./tests/blast-perf/run.sh cores            # single invocation, prints run_id
./tests/blast-perf/run.sh cores-report <run_id>
```

Everything runs in one `cores.nf` invocation now (cold pass + full warm
sweep on one node), so there's no scale-down wait and no second mode —
unlike the `batch`/`per_query` split, there's only one shape being varied
here (thread count), not two things that need isolated nodes.

## Report: `summarise_cores.py` → `REPORT.md`

A new script rather than more branches in `summarise.py`. It reuses the
helpers from there by importing them (`read_trace`, `read_timing_files`,
`parse_nextflow_duration`, `format_duration`). Standard library only, run
from `./venv`, flake8 clean.

| Metric | cold (48) | warm 2 | warm 4 | warm 8 | warm 16 | warm 48 |
|---|---|---|---|---|---|---|
| blastn time | | | | | | |
| CPU-seconds (from samples, not just trace avg) | | | | | | |
| Peak concurrent threads (from samples) | | | | | | |
| `cpu.max`, `memory.max` | (once, shared across all passes on this node) |

Plus workflow walltime, NODE_UP duration + hostname, and estimated cost,
as in the first report.

**From the samples:** compute instantaneous CPU-seconds/sec between
consecutive sample lines, so the report can say whether a pass actually
used N cores at any point, or just the average implies it.

**Checks:**

- **Correctness:** every pass's `.bls` (cold and all five warm) sorted and
  diffed against each other. Thread count shouldn't change MEGABLAST
  hits. Report any differing row count.
- **Warm really was warm:** flag it if any warm pass's disk-sector delta
  (from the samples) is anywhere close to the cold pass's — it would mean
  cache eviction happened mid-sweep, invalidating later points.

**Verdict:** a table of blastn time vs thread count (warm), stating the
speedup of each step relative to 2 threads, plus the cold-pass number on
its own with a note that only one cold data point exists. Also give
**throughput per node**: how many BLAST tasks at each thread count could
run concurrently on a 48-core/384GB node, times that pass's speed, to show
which thread count gets the most total work done per node-hour. State
plainly that this is one run, so small differences between adjacent
thread counts are inconclusive.

## Implementation order

1. `BLAST_PASS` (with sampling loop) in `modules.nf`, `cores.nf`, config
   blocks.
2. Smoke-test locally against the test DB — confirms wiring and the
   sampling loop's fallback path (no real NVMe device locally).
3. `summarise_cores.py` against the smoke outputs; flake8.
4. `run.sh` mode; README section.
5. Azure: `run.sh cores`, then `run.sh cores-report <run_id>`.

## Cost and time estimate

One invocation: ~15–20 minutes node start-up + refdata staging, plus the
cold pass (~2–3 minutes, per the last run) and five warm passes. Warm
passes at 611 CPU-seconds of work (the batch total, not per-query):
2 threads ≈ 5min, 4 ≈ 2.5min, 8 ≈ 1.3min, 16 ≈ 40s, 48 ≈ 15–30s (thread
overhead likely keeps it off the theoretical floor at the high end). Sum
of warm passes is roughly 10–12 minutes. Plus ~15 minutes keep-warm at the
end.

Total: ~45–55 minutes on one node, roughly **$6–7**. Cheaper than the
original 4-vs-48 two-invocation plan ($11–13), since this needs only one
node instead of two.

## Follow-ups (outside this test)

1. If low thread counts win on throughput per node, test the shape that
   follows: many per-query (or small-batch) tasks at that thread count,
   all running at once on one node. That's the production shape this
   result would point to.
2. If `cpu.max` shows no limit, the production `cpus` setting is only a
   scheduling hint, and `-num_threads` (driven by `blast_threads` in
   `params/default_params.yml`, currently 2) is what actually controls
   BLAST (see follow-up 1 in `batch-vs-series.md`).
