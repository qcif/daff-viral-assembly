#!/usr/bin/env python
import pandas as pd
import argparse
import os
from functools import reduce
import pytaxonkit
import numpy as np
import re
from pytaxonkit import lineage


def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Load taxonomic classification and summarise.")
    parser.add_argument("--kaiju", required=True, type=str)
    parser.add_argument("--sample_name", required=True, type=str)
    parser.add_argument("--bracken", required=True, type=str)
    parser.add_argument("--taxonkit_database_dir", required=True, type=str)
    parser.add_argument("--stats", required=True, type=str, help="path to bbsplit stats file")
    parser.add_argument("--filter", required=True, type=str, help="path to file with terms to filter out")
    return parser.parse_args()

def load_blast_results(path):
    """Load BLASTn results based on mode."""
    if not os.path.isfile(path):
        raise FileNotFoundError(f"BLASTn results file not found: {path}")

    # Define expected header (based on your BLAST output fields)
    columns = [
        "qseqid", "sgi", "sacc", "alignment_length", "pident", "mismatch", "gapopen",
        "qstart", "qend", "qlen", "sstart", "send", "slen", "sstrand", "evalue",
        "bitscore", "qcovhsp", "stitle", "staxids", "qseq", "sseq", "sseqid",
        "qcovs", "qframe", "sframe"
    ]

    df = pd.read_csv(path, sep="\t", header=None, dtype=str)
    if df.shape[1] != len(columns):
        raise ValueError(
            f"Expected {len(columns)} columns (BLAST output), but found {df.shape[1]} in {path}"
        )
    df.columns = columns

    # Define columns and dtypes
    dtype = {
        "qseqid": 'str', "sgi": 'str', "sacc": 'str', "alignment_length": 'int64', "pident": 'float64', "mismatch": 'int64',
        "gapopen": 'int64', "qstart": 'int64', "qend": 'int64', "qlen": 'int64', "sstart": 'int64',
        "send": 'int64', "slen": 'int64', "sstrand": 'str', "evalue": 'float64', "bitscore": 'float64',
        "qcovhsp": 'int64', "stitle": 'str', "staxids": 'str', "qseq": 'str', "sseq": 'str', "sseqid": 'str',
        "qcovs": 'int64', "qframe": 'int64', "sframe": 'int64'
    }
    # Load DataFrame
    for col, dtype_ in dtype.items():
        if col in df.columns:
            if dtype_ in ("int64", "float64"):
                df[col] = pd.to_numeric(df[col], errors="coerce").astype(dtype_)
            else:
                df[col] = df[col].astype(str)

    df["staxids"] = pd.to_numeric(df["staxids"].str.split(";").str[0], errors='coerce').fillna(0).astype(int)
    return df

def enrich_with_taxonomy(df, taxonkit_dir):
    """Add taxonomy information to the DataFrame."""
    staxids_l = df["staxids"].unique().tolist()
    lineage_df = pytaxonkit.lineage(staxids_l, data_dir=taxonkit_dir)[['TaxID', 'FullLineage']]
    lineage_df.columns = ["staxids", "full_lineage"]
    lineage_df["staxids"] = lineage_df["staxids"].astype(int)
    lineage_df["full_lineage"] = lineage_df["full_lineage"].str.lower().str.replace(" ", "_", regex=False)
    # ensure no NaNs and treat values as strings
    lineage_df["full_lineage"] = lineage_df["full_lineage"].fillna("NA").astype(str)
    lineage_df["broad_taxonomic_category"] = np.where(
        lineage_df["full_lineage"].str.contains("viruses;"),
        "virus",
        np.where(
            lineage_df["full_lineage"].str.contains("viroids;"),
            "viroids",
            "non-viral"
        )
    )

    names_df = pytaxonkit.name(staxids_l, data_dir=taxonkit_dir)[['TaxID', 'Name']]
    names_df.columns = ["staxids", "species"]
    names_df["staxids"] = names_df["staxids"].astype(int)

    return [df, names_df, lineage_df]

def merge_taxonomy(dfs):
    """Merge taxonomy-enriched data."""
    return reduce(lambda left, right: pd.merge(left, right, on="staxids", how="outer"), dfs)

def filter_and_format(df, sample_name, filter_file):
    """Filter only viral hits and format."""
    df.insert(0, "sample_name", sample_name)
    df = df[~df["species"].str.contains("synthetic construct", na=False)]
    df = df[~df["species"].str.contains("Expression vector", na=False)]
    # Then drop duplicates — keep first valid hit per qseqid
    df = df.sort_values(by=["qseqid", "bitscore"], ascending=[True, False])  # Optional: Sort by best score
    df = df.drop_duplicates(subset=["qseqid"], keep="first").copy()


    # Filter only viral entries
    df = df[
        df["broad_taxonomic_category"].notna()
        & (df["broad_taxonomic_category"].str.strip() != "")
    ].copy()
    df = df[df["broad_taxonomic_category"] == "virus"].copy()


    df["RNA_type"] = np.where(
        df.stitle.str.contains("RNA1|RNA 1|segment 1|polyprotein P1", case=False, na=False), "RNA1",
        np.where(
            df.stitle.str.contains("RNA2|RNA 2|segment 2|polyprotein P2", case=False, na=False), "RNA2",
            np.where(
                df.stitle.str.contains("RNA3|RNA 3|segment 3|polyprotein P3", case=False, na=False), "RNA3",
                ""
            )
        )
    )

    df["species_updated"] = df[["species", "RNA_type"]].fillna("").astype(str).agg(" ".join, axis=1)
    df["species_updated"] = df["species_updated"].str.replace("RNA1 RNA1", "RNA1", regex=False)
    df["species_updated"] = df["species_updated"].str.replace("RNA2 RNA2", "RNA2", regex=False)
    df["species_updated"] = df["species_updated"].str.replace("RNA3 RNA3", "RNA3", regex=False)
    df["species_updated"] = df["species_updated"].str.rstrip()

    df["ncontigs"] = df.groupby(["species_updated", "sacc"])["sacc"].transform("count")
    df["ncontigs_score"] = df.groupby("species_updated")["ncontigs"].transform(max_naccs)

    df["pident"] = df["pident"].astype(float)
    df["pident_score"] = df.groupby("species_updated")["pident"].transform(max_pid)

    df["bitscore"] = df["bitscore"].astype(int)
    df["bitscore_score"] = df.groupby("species_updated")["bitscore"].transform(max_bitscore)

    df["evalue"] = df["evalue"].astype(float)
    df["evalue_score"] = df.groupby("species_updated")["evalue"].transform(min_evalue)

    # Extract cov value using regex and convert to float
    df["assembly_kmer_cov"] = df["qseqid"].str.extract(r'_cov_([0-9.]+)').astype(float)
    df["assembly_kmer_cov_score"] = df.groupby("species_updated")["assembly_kmer_cov"].transform(max_assembly_kmer_cov)

    # assign a score of 1 if qcovs > 75, else 0
    df["qcovs_score"] = global_qcovs(df)
    df["best_qcovs_score"] = df.groupby("species_updated")["qcovs"].transform(max_qcovs)
    df["completeness_score"] = df["stitle"].apply(completeness_score)

    df["qlen"] = df["alignment_length"].astype(int)
    df["query_length_score"] = df.groupby("species_updated")["qlen"].transform(max_length)
    df["alignment_length"] = df["alignment_length"].astype(int)
    df["alignment_length_score"] = df.groupby("species_updated")["alignment_length"].transform(max_length)

    df["total_score"] = (
        df["ncontigs_score"]
        + df["pident_score"]
        + df["bitscore_score"]
        + df["evalue_score"]
        + df["assembly_kmer_cov_score"]
        + df["qcovs_score"]
        + df["best_qcovs_score"]
        + df["completeness_score"]
        + df["alignment_length_score"]
        + df["query_length_score"]
    )

    # Find the index of the best hit per species
    best_idx = df.groupby("species_updated")["total_score"].idxmax()

    #Create a new column indicating whether the row is the best hit
    df["best_contig_per_sp_filter"] = df.index.isin(best_idx)


    # Read exclusion patterns from a file
    with open(filter_file, "r") as f:
        exclude_patterns = [line.strip() for line in f if line.strip()]

    pattern = "|".join(exclude_patterns)

    # Create a new column that indicates whether the row is filtered
    df["term_filter"] = ~(
        df["species_updated"].str.contains(pattern, case=False, na=False)
        | df["stitle"].str.contains(pattern, case=False, na=False)
        )

    df["cov_filter"] = (
        (df["assembly_kmer_cov"] >= 5) &
        (df["qcovs"] >= 30)
        )



    final_columns_filt = ["sample_name", "qseqid", "sacc", "alignment_length", "evalue", "bitscore", "pident", "mismatch",
                     "gapopen", "qstart", "qend", "qlen", "sstart", "send", "slen", "sstrand",
                     "qcovhsp", "staxids", "qseq", "sseq", "qcovs",
                     "species_updated", "RNA_type", "stitle", "full_lineage", "ncontigs", "total_score", "term_filter", "cov_filter", "best_contig_per_sp_filter"]

    return df[final_columns_filt]

def parse_bbsplit_log(logfile):
    total_reads = None
    r1_mapped = r2_mapped = 0

    with open(logfile, "r") as f:
        section = None
        for line in f:
            line = line.strip()

            # total input reads
            m_total = re.match(r"Reads Used:\s+(\d+)", line)
            if m_total:
                total_reads = int(m_total.group(1))
                continue

            # detect sections
            if line.startswith("Read 1 data"):
                section = "R1"
                continue
            elif line.startswith("Read 2 data"):
                section = "R2"
                continue

            # mapped reads per read
            m_mapped = re.match(r"mapped:\s+[\d\.%]+\s+(\d+)", line)
            if m_mapped:
                count = int(m_mapped.group(1))
                if section == "R1":
                    r1_mapped = count
                elif section == "R2":
                    r2_mapped = count

    if total_reads is None:
        raise ValueError("Could not find total reads in log file.")

    total_mapped_reads = r1_mapped + r2_mapped
    clean_reads = total_reads - total_mapped_reads



    return clean_reads

def is_viral(tid, taxid_to_lineage):
    if not tid or tid == 0:
        return False
    lineage_path = taxid_to_lineage.get(tid, "")
    if not lineage_path:
        return False
    lineage_ids = [x for x in lineage_path.split(";") if x]
    return ("10239" in lineage_ids) or ("12884" in lineage_ids)

def categorize(row, taxid_to_lineage):
    name = str(row.get("taxon_name", "") or "").lower()
    lineage = str(row.get("full_lineage", "") or "").lower()

    # Generalized unclassified virus check
    if "unclassified" in name and ("virus" in name or "virinae" in name):
        return "unclassified virus"

    # Split lineage into parts
    lineage_parts = [x.strip() for x in lineage.split(";") if x.strip()]
    last_part = lineage_parts[-1] if lineage_parts else ""

    # Check if lineage contains 'viruses' to identify virus clade
    is_viral_lineage = "viruses" in lineage_parts[0] if lineage_parts else False

    # Categorize
    if is_viral_lineage and last_part.startswith("unclassified"):
        return "unclassified virus"
    elif is_viral_lineage:
        return "viral"
    elif "unclassified" in name:
        return "unclassified"
    elif "cannot be assigned to a (non-viral) species" in name:
        return "unclassified virus"
    else:
        return "other non-viral"

def categorize_bracken(row, taxid_to_lineage):
    tid = row["taxon_id"]
    if is_viral(tid, taxid_to_lineage):
        return "viral"
    else:
        return "other non-viral"
def chunked(iterable, size=500):
    """Yield successive chunks from a list."""
    for i in range(0, len(iterable), size):
        yield iterable[i:i + size]

def main():
    args = parse_arguments()
    kaiju_path = args.kaiju
    sample_name = args.sample_name
    bracken_path = args.bracken
    tk_db_dir = args.taxonkit_database_dir
    log = args.stats
    filter_file = args.filter


    with open(filter_file, "r") as f:
        exclude_patterns = [line.strip() for line in f if line.strip()]

    pattern = "|".join(exclude_patterns)

    if not os.path.isfile(kaiju_path):
        raise FileNotFoundError(f"{kaiju_path} does not exist.")

    if not os.path.isfile(bracken_path):
        raise FileNotFoundError(f"{bracken_path} does not exist.")
    df = pd.read_csv(kaiju_path, sep="\t", dtype=str)


    # ------------------------------------
    # Load Kaiju summary table
    # ------------------------------------

    # Convert taxon_id safely
    # Get lineage only for non-zero taxids
    unique_taxids = df.loc[df["taxon_id"] != 0, "taxon_id"].unique().tolist()
    print(f"Unique taxids in Kaiju results: {len(unique_taxids)}")

    # Query lineage (just pass taxid list)
    lin = lineage(unique_taxids, data_dir=tk_db_dir)
    
    # Convert to DataFrame
    lin_df = pd.DataFrame(lin)
    print(lin_df.columns)
    print(lin_df.iloc[0])
    # Keep relevant columns and rename
    lin_df = lin_df.rename(columns={
        "TaxID": "taxid",
        "FullLineageTaxIDs": "full_lineage_taxids",
        "FullLineage": "full_lineage_names",
        "FullLineageRanks": "full_lineage_ranks",
        "Name": "name"
    })
    lin_df["full_lineage_names"] = lin_df["full_lineage_names"].fillna("NA").astype(str)
    lin_df["full_lineage_ranks"] = lin_df["full_lineage_ranks"].fillna("NA").astype(str)
    # Fill NaNs with 0 and convert to int
    lin_df["taxid"] = pd.to_numeric(lin_df["taxid"], errors="coerce").fillna(0).astype(int)



    # Build mapping
    taxid_to_lineage = {int(k): str(v) for k, v in zip(lin_df["taxid"], lin_df["full_lineage_taxids"])}

    taxid_to_lineage_names = dict(zip(
        lin_df["taxid"],
        lin_df["full_lineage_names"]
    ))
    taxid_to_lineage_ranks = dict(zip(
        lin_df["taxid"],
        lin_df["full_lineage_ranks"]
    ))

    # Viral check
    df["taxon_id"] = pd.to_numeric(df["taxon_id"], errors="coerce").fillna(0).astype(int)

    df["full_lineage"] = df["taxon_id"].apply(lambda x: taxid_to_lineage_names.get(int(x), ""))
    df["full_lineage_ranks"] = df["taxon_id"].apply(
        lambda x: taxid_to_lineage_ranks.get(int(x), "")
    )
    df["broad_categories"] = df.apply(
        lambda row: categorize(row, taxid_to_lineage),
        axis=1
    )
    df.loc[df["taxon_name"] == "cannot be assigned to a (non-viral) species",
       ["full_lineage", "full_lineage_ranks"]] = "NA"
    df.loc[df["taxon_name"] == "unclassified",
       ["full_lineage", "full_lineage_ranks"]] = "NA"

    df["term_filter"] = ~(
        df["taxon_name"].str.contains(pattern, case=False, na=False)
        )

    df = df.rename(columns={
        "percent": "pc_reads"
    })
    df["cov_filter"] = df["pc_reads"].astype(float) >= 0.002
    df_subset = df[['taxon_name', 'taxon_id', 'full_lineage', 'full_lineage_ranks', 'broad_categories', 'reads', 'pc_reads', 'term_filter', 'cov_filter']]
    df_subset.to_csv(sample_name + "_kaiju_summary.txt", sep="\t", index=False)

    # ------------------------------------
    # Load KRAKEN/BRAKEN summary table
    # ------------------------------------
    br = pd.read_csv(bracken_path, sep="\t", dtype=str)
    br = br.rename(columns={
        "taxonomy_id": "taxon_id"})
    unique_taxids_br = br.loc[br["taxon_id"] != 0, "taxon_id"].unique().tolist()
    lin2 = lineage(unique_taxids_br, data_dir=tk_db_dir)
    lin2_df = pd.DataFrame(lin2)
    lin2_df = lin2_df.rename(columns={
        "TaxID": "taxid",
        "FullLineageTaxIDs": "full_lineage_taxids",
        "FullLineage": "full_lineage_names",
        "FullLineageRanks": "full_lineage_ranks",
        "Name": "name"
    })
    lin2_df["full_lineage_ranks"] = lin2_df["full_lineage_ranks"].fillna("").astype(str)
    lin2_df["taxid"] = pd.to_numeric(lin2_df["taxid"], errors="coerce").fillna(0).astype(int)
    lin2_df["full_lineage_taxids"] = lin2_df["full_lineage_taxids"].fillna("")
    taxid_to_lineage = {int(k): str(v) for k, v in zip(lin2_df["taxid"], lin2_df["full_lineage_taxids"])}
    taxid_to_lineage_names = dict(zip(
        lin2_df["taxid"],
        lin2_df["full_lineage_names"]
    ))
    taxid_to_lineage_ranks = dict(zip(
        lin2_df["taxid"],
        lin2_df["full_lineage_ranks"]
    ))
    br["taxon_id"] = pd.to_numeric(br["taxon_id"], errors="coerce").fillna(0).astype(int)
    br["full_lineage"] = br["taxon_id"].apply(lambda x: taxid_to_lineage_names.get(int(x), ""))
    br["full_lineage_ranks"] = br["taxon_id"].apply(
        lambda x: taxid_to_lineage_ranks.get(int(x), "")
    )
    br["broad_categories"] = br.apply(
        lambda row: categorize_bracken(row, taxid_to_lineage),
        axis=1
    )

    filtered_read_counts = parse_bbsplit_log(log)
    br["new_est_reads"] = pd.to_numeric(br["new_est_reads"], errors="coerce")
    br = br.rename(columns={
        "name": "taxon_name",
        "new_est_reads": "reads",
    })
    br["term_filter"] = ~(
        br["taxon_name"].str.contains(pattern, case=False, na=False)
        )
    # Sum of all classified reads in Bracken table
    classified_reads = br["reads"].sum()
    # Compute unclassified reads
    unclassified_reads = filtered_read_counts - classified_reads
    unclassified_reads = max(unclassified_reads, 0)  # safety
    # Build the unclassified row
    unclassified_row = {
        "taxon_name": "unclassified",
        "taxon_id": 0,
        "full_lineage": "NA",
        "full_lineage_ranks": "NA",
        "broad_categories": "unclassified",
        "reads": unclassified_reads,
        "term_filter": True
    }
    # Append to the bottom of the table
    br = pd.concat([br, pd.DataFrame([unclassified_row])], ignore_index=True)
    br["pc_reads"] = (br["reads"] / filtered_read_counts) * 100
    br["cov_filter"] = br["pc_reads"].astype(float) >= 0.001
    br_subset = br[['taxon_name', 'taxon_id', 'full_lineage', 'full_lineage_ranks', 'broad_categories', 'reads', 'pc_reads', 'term_filter', 'cov_filter']]
    br_subset.to_csv(sample_name + "_kraken_summary.txt", sep="\t", index=False)

if __name__ == "__main__":
    main()

#singularity exec  -B /work/daff_viral_rnaseq/rnaspades/BR_Guar_18 docker://quay.io/biocontainers/pytaxonkit:0.9.1--pyhdfd78af_1 /work/daff_viral_rnaseq/rnaspades/BR_Guar_18/filter_blast.py --blastn_results BR_Guar_18_blastn_header.txt --sample_name BR_Guar_18 --mode ncbi --taxonkit_database_dir ~/.taxonkit



#--sample_name Dataset_1 --bracken /work/Dataset_1_bracken_report.txt --kaiju /work/Dataset_1_kaiju_summary.tsv --taxonkit_database_dir ~/.taxonkit --stats /work/Dataset_1_bbsplit_stats.txt --filter /work/bin/filterKeyWords.txt