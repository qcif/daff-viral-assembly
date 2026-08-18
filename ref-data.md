# Reference data for VIEW

## Provenance

The canonical names discard version information, so they recorded it here.

| Path | Source | Version / build |
|---|---|---|
| `kraken2_db/` | `k2_core_nt_20251015.tar.gz` | core_nt, 2025-10-15 |
| `kaiju_db/kaiju_db.fmi` | `kaiju_db_nr_euk_2023-05-10.tgz` | nr_euk, 2023-05-10 |
| `pfam/Pfam-A.hmm` | Pfam `current_release` | as of 2026-07-21 |
| `diamond/viral.dmnd` | `U-RVDBv32.0-prot.fasta.xz` | RVDB v32.0 (protein) |
| `rvdb_taxonomy` | `RVDB_Taxon_Current.tab.gz` | as of 2026-07-21 |
| `genomad_db/` | Zenodo record 14886553 | v1.9 |
| `rrna_ref` | SortMeRNA `database.tar.gz` | v4.3.4 sensitive, U→T converted |


## Build the reference data

```sh
# This worked as of 2026-07-21, but will likely need to be updated in future

# Build rvdb.dmnd
wget https://rvdb-prot.pasteur.fr/files/U-RVDBv32.0-prot.fasta.xz
unxz U-RVDBv32.0-prot.fasta.xz
diamond makedb --quiet --threads 2 --in U-RVDBv32.0-prot.fasta -d rvdb
rm U-RVDBv32.0-prot.fasta

# RVDB taxonomy
wget https://rvdb.dbi.udel.edu/download/RVDB_Taxon_Current.tab.gz

# Kraken
wget https://genome-idx.s3.amazonaws.com/kraken/k2_core_nt_20251015.tar.gz

# Kaiju
wget https://kaiju-idx.s3.eu-central-1.amazonaws.com/2023/kaiju_db_nr_euk_2023-05-10.tgz

# GeNomad DB
wget -O genomad_db_v1.9.tar.gz 'https://zenodo.org/records/14886553/files/genomad_db_v1.9.tar.gz?download=1'

# HMMER PFAM-A.hmm
# NOTE: Pfam-A.hmm.gz is the HMM library. Pfam-A.hmm.dat.gz is only the
# annotation flatfile and cannot be used by hmmpress/hmmscan.
wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz

# SortMeRNA rRNA references
wget -O sortmerna_database.tar.gz https://github.com/sortmerna/sortmerna/releases/download/v4.3.4/database.tar.gz
```

Extract and decompress:

```sh
tar -xzf k2_core_nt_20251015.tar.gz -C k2_core_nt
tar -xzf kaiju_db_nr_euk_2023-05-10.tgz -C kaiju
tar -xzf genomad_db_v1.9.tar.gz
tar -xzf sortmerna_database.tar.gz
gunzip RVDB_Taxon_Current.tab.gz
gunzip Pfam-A.hmm.gz
```

The Kaiju tarball ships `names.dmp` and `nodes.dmp` alongside the `.fmi`. Keep
them — the pipeline globs the DB directory for all three, and the taxonomy must
match the index build date (2023-05-10). Do not substitute the `.dmp` files from
the Kraken2 or geNomad databases; mismatched taxon IDs misclassify silently.

Index the Pfam HMM library. All four output files must sit beside the `.hmm` and
share its basename:

Via the same container the pipeline uses:

```sh
docker run --rm -u "$(id -u):$(id -g)" -v "$PWD:/db" \
    quay.io/biocontainers/hmmer:3.4--hb6cb901_4 \
    hmmpress /db/Pfam-A.hmm
```

Convert the SortMeRNA reference from RNA to DNA. BBDuk matches DNA reads, so
uracils in the reference silently prevent rRNA removal:

```sh
sed '/^>/! y/Uu/Tt/' smr_v4.3_sensitive_db.fasta > smr_v4.3_sensitive_db_DNA.fasta
```

Verify the conversion — sequence lines must contain no uracil. The check ignores
headers, which legitimately contain "U":

```sh
grep -v '^>' smr_v4.3_sensitive_db_DNA.fasta | grep -c '[Uu]'   # must print 0
```

## Final layout on disk

Rename to the canonical names. This layout is shared by the data server, the
`refdata-wf4` blob container, and `/mnt/nvme/refdata/` on the Batch nodes, and
matches the paths hardcoded in `deploy/azure/run-wf4.sh`:

```sh
mv k2_core_nt kraken2_db
mv kaiju kaiju_db
mv kaiju_db/kaiju_db_nr_euk.fmi kaiju_db/kaiju_db.fmi

mkdir -p pfam diamond
mv Pfam-A.hmm* pfam/
mv rvdb.dmnd diamond/viral.dmnd

mv RVDB_Taxon_Current.tab rvdb_taxonomy
mv smr_v4.3_sensitive_db_DNA.fasta rrna_ref
```

Only the converted rRNA reference is retained. The remaining SortMeRNA fastas
(`default`, `fast`, `rfam_seeds`, and the unconverted `sensitive`) are unused by
the pipeline and are deliberately not kept, so that what is on disk matches what
is staged to the nodes.

Resulting tree:

```
.
├── diamond
│   └── viral.dmnd
├── genomad_db
│   ├── genomad_db
│   ├── genomad_db.dbtype
│   ├── genomad_db_h
│   ├── genomad_db_h.dbtype
│   ├── genomad_db_h.index
│   ├── genomad_db.index
│   ├── genomad_db.lookup
│   ├── genomad_db_mapping
│   ├── genomad_db.source
│   ├── genomad_db_taxonomy
│   ├── genomad_integrase_db
│   ├── genomad_integrase_db.dbtype
│   ├── genomad_integrase_db_h
│   ├── genomad_integrase_db_h.dbtype
│   ├── genomad_integrase_db_h.index
│   ├── genomad_integrase_db.index
│   ├── genomad_integrase_db.lookup
│   ├── genomad_integrase_db.source
│   ├── genomad_marker_metadata.tsv
│   ├── genomad_mini_db -> genomad_db
│   ├── genomad_mini_db.dbtype
│   ├── genomad_mini_db_h -> genomad_db_h
│   ├── genomad_mini_db_h.dbtype -> genomad_db_h.dbtype
│   ├── genomad_mini_db_h.index -> genomad_db_h.index
│   ├── genomad_mini_db.index
│   ├── genomad_mini_db.lookup -> genomad_db.lookup
│   ├── genomad_mini_db_mapping -> genomad_db_mapping
│   ├── genomad_mini_db.source -> genomad_db.source
│   ├── genomad_mini_db_taxonomy -> genomad_db_taxonomy
│   ├── mini_set_ids
│   ├── names.dmp
│   ├── nodes.dmp
│   ├── plasmid_hallmark_annotation.txt
│   ├── version.txt
│   └── virus_hallmark_annotation.txt
├── kaiju_db
│   ├── kaiju_db.fmi
│   ├── names.dmp
│   └── nodes.dmp
├── kraken2_db
│   ├── database100mers.kmer_distrib
│   ├── database150mers.kmer_distrib
│   ├── database200mers.kmer_distrib
│   ├── database250mers.kmer_distrib
│   ├── database300mers.kmer_distrib
│   ├── database50mers.kmer_distrib
│   ├── database75mers.kmer_distrib
│   ├── hash.k2d
│   ├── inspect.txt
│   ├── ktaxonomy.tsv
│   ├── library_report.tsv
│   ├── names.dmp
│   ├── nodes.dmp
│   ├── opts.k2d
│   ├── seqid2taxid.map
│   └── taxo.k2d
├── pfam
│   ├── Pfam-A.hmm
│   ├── Pfam-A.hmm.h3f
│   ├── Pfam-A.hmm.h3i
│   ├── Pfam-A.hmm.h3m
│   └── Pfam-A.hmm.h3p
├── rrna_ref
└── rvdb_taxonomy
```
