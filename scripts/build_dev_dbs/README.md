# Dev reference database builders

Scripts to build small (~500 MB each) local Kraken2 and Kaiju databases for
development iteration on the wf4 pipeline. These are **not** the CI smoke-test
DBs under [`tests/`](../../tests/) (those are ~40 KB and classify essentially
nothing), and they are **not** the full production DBs (~490 GB combined).

## Contents

| File | Purpose |
|---|---|
| `seed_taxa.txt` | Shared seed list of taxa (viruses / fabaceae / context), grouped into `[section]` blocks. Edit this to change what both DBs cover. |
| `build_kraken2_dev_db.sh` | Builds `refdata_dev/kraken2_db/` from `seed_taxa.txt` + viral RefSeq. |
| `build_kaiju_dev_db.sh` | Builds `refdata_dev/kaiju_db/kaiju_dev_db.fmi` from `seed_taxa.txt` + viral RefSeq protein. |
| `lib.sh` | Shared bash helpers (seed parsing, NCBI E-utilities fetches). |

Both scripts write outputs relative to the repo root (`refdata_dev/`, gitignored).

## Prerequisites

- `docker` — kraken2 and kaiju are pulled as containers; nothing else to install
- `curl`, `awk`, `tar`, `gzip` (standard)

Default images (override via env if needed):

- `staphb/kraken2:2.1.3` — set `KRAKEN2_IMAGE` to change
- `nanozoo/kaiju:1.10.1--90efeeb` — set `KAIJU_IMAGE` to change

## Usage

From the repo root (writes to `refdata_dev/`):

```sh
bash scripts/build_dev_dbs/build_kraken2_dev_db.sh
bash scripts/build_dev_dbs/build_kaiju_dev_db.sh
```

Or to write elsewhere (e.g. a large scratch disk):

```sh
OUT_ROOT=/mnt/data_2/dev_dbs bash scripts/build_dev_dbs/build_kraken2_dev_db.sh
OUT_ROOT=/mnt/data_2/dev_dbs bash scripts/build_dev_dbs/build_kaiju_dev_db.sh
```

Both scripts are idempotent — rerun to resume after a failure. Delete the
target directory (`refdata_dev/kraken2_db` or `refdata_dev/kaiju_db`) to force
a full rebuild.

## Expected cost

|  | Runtime | Peak working disk | Final size |
|---|---|---|---|
| Kraken2 | ~15–25 min | ~5 GB | ~500 MB |
| Kaiju | ~10–15 min | ~3 GB | ~500 MB |

Runtime is dominated by NCBI downloads on the first run; subsequent runs use
cached files under `refdata_dev/.build/`.

## Pointing the workflow at the dev DBs

A `dev` Nextflow profile is wired up in [`nextflow.config`](../../nextflow.config)
to point `params.kraken2_db` and `params.kaiju_db` at `refdata_dev/`. Run with:

```sh
nextflow run main.nf -profile singularity,dev \
  --input tests/index_test.csv \
  --analyst_name "Dev" --facility "Local"
```

## Editing the seed list

`seed_taxa.txt` is grouped into three sections:

- `[viruses]` — currently the whole viral RefSeq (taxid 10239). Rarely needs
  editing.
- `[fabaceae]` — plant hosts. Must include *Clitoria ternatea* (`NC_042239.1`).
  Add legume relatives as new accessions on their own line.
- `[context]` — generic non-target genomes so host/contaminant reads have
  somewhere to land.

Entries are either numeric NCBI taxids or NCBI nucleotide accessions
(`NC_...`, `NZ_...`, etc.). Blank lines and `#` comments are ignored.
