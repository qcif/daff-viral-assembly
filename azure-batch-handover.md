# VIEW on Azure Batch — test setup

Branch `azure-refdata-test`.

## The pool

`view_test` — a single `Standard_L8as_v3` (8 cores, 8 task slots), autoscaling
0→1 on pending tasks. Its start task
(`deploy/azure/setup-wf4-test.sh.template`, uploaded to
`daffstandard/scripts/`) stages the reference data that can't be shipped from
the launcher: `taxdump` and `genomad_db`, both to `/mnt/nvme/refdata/`. Everything
else the test case needs is the ~8MB of test DBs in `tests/`, which Nextflow
stages per-task as normal `path` inputs.

Two notes on that start task:

- The taxdump is **doubly nested** on blob (`refdata/taxdump/taxdump/*.dmp`) — the
  inner dir is the real source.
- geNomad's `genomad_mini_db_*` symlinks can't live on blob, so they're excluded
  from the upload and recreated on-node. `genomad_mini_db.dbtype` and `.index` are
  *not* symlinks upstream — they're real files and come down with the rest.

Batch allows 2 pools and `taxodactyl` holds one, so `view` was deleted to make
room. Prod and test can't currently coexist.

## Entrypoint

```sh
./deploy/azure/run-wf4-test.sh [-resume]
```

Runs `-profile azure_test,test`. The trailing `test` supplies the `tests/` DB
paths, so it has to stay last.

## Code changes

**`path` → `val` for node-staged data.** Nextflow resolves `path` inputs on the
*launching* machine and stages them; there's no "already exists on the executor"
option. So anything staged by the start task has to be a `val`, passed as a plain
string, with the location bind-mounted into the container
(`containerOptions` in `conf/azure.test.config`). Converted: `taxdump` in
`summarise_read_classification`, `extract_raw_viral_blast_hits`,
`extract_final_viral_blast_hits`; `genomad_db` in `genomad/endtoend` plus
`Channel.value` in `workflows/view/main.nf`. The local docker profile
bind-mounts the same paths via a lazy `containerOptions` closure.

**geNomad no longer self-downloads.** `GENOMAD_DOWNLOAD_DB` fires when
`genomad_db` is null and the profile contains `test` — it pulled ~750MB per run
and published it back to the launcher via `publishDir "${params.databases}"`,
leaving an un-gitignored 733MB dir in the repo. `conf/azure.test.config` now sets
`genomad_db` to the node path so the branch never fires on Azure; `databases/` is
gitignored. The download branch still works locally.

**Container images need a registry host.** Batch rejects unqualified names
(`One or more container images specified are invalid`). Six were wrong —
`CAT_FASTQ`, `FASTQC_TRIM`, `FQ_SUBSAMPLE`, `KAIJU_KAIJU`, `SEQTK_SAMPLE` used
Docker Hub shorthand; `HTML_REPORT` used a `docker://` URI. Fixed in
`nextflow.config`, re-pointed to `quay.io` at identical tags.

**`withName` selectors.** They match the *aliased* name from
`workflows/view/main.nf` as a full-string regex, and silently no-op otherwise.
Fixed `KRAKEN2`→`KRAKEN2_KRAKEN2`, `KAIJU`→`KAIJU_KAIJU`, `BWA`→`MAP_TO_CONTIGS`,
`CDHIT_CDHIT`→`CLUSTER`, `SEQTK_SUBSEQ`→`EXTRACT_CONTIGS`. The one that mattered
was `maxForks = 1` on Kraken2/Kaiju — unapplied, two concurrent 300GB-label tasks
OOM the node.

**Watch out:** identical `withName` keys across config files **replace** rather
than merge, silently dropping `container`. That's caused three outages so far
(`MEGABLAST*`, `HMMSCAN`, `KAIJU_KAIJU`); colliding blocks now restate their
container. Any future resource-only override in an Azure config will re-trigger it.

**`errorStrategy`** `ignore` → `terminate`, so failures stop the run instead of
surfacing as missing outputs.

**Autoscale formula** (`deploy/azure/batch-helpers.sh`): `GetSample()` aborts the
whole formula when a metric has no history, and Batch falls back to 0 nodes —
deadlocking a fresh pool. Every call is now guarded by `GetSamplePercent()`.

## Not working yet

The test case reaches `SUMMARISE_READ_CLASSIFICATION` (16 processes pass) and
fails there. Six params default to `${projectDir}`/`${launchDir}` paths and are
interpolated raw into process scripts, so they only exist on the launcher:
`filter_terms` (4 modules), `tool_versions`, `default_params`, `phix`. These are
small files, so the fix is the inverse of the above — make them real `path`
inputs and let Nextflow stage them.

`BBMAP_BBSPLIT` passed despite `params.phix` not existing. If `bbsplit.sh`
tolerates a missing `ref=`, phiX removal was silently skipped — worth checking
`.command.err` before trusting that output.
