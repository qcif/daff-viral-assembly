#!/usr/bin/env python

import pandas as pd
import argparse
import os
from functools import reduce
import pytaxonkit
import numpy as np
import tempfile
import shutil
import re

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Load and enrich BLASTn results and select only viruses.")
    parser.add_argument("--blastn_results", required=True, type=str)
    parser.add_argument("--sample_name", required=True, type=str)
    parser.add_argument("--taxonkit_database_dir", required=True, type=str)
    parser.add_argument("--filter", required=True, type=str)
    parser.add_argument("--assembly_headers", required=True, type=str)
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
    # Create a temporary file with header + original content
    df = pd.read_csv(path, sep="\t", header=None, dtype=str)
    if df.shape[1] != len(columns):
        raise ValueError(
            f"Expected {len(columns)} columns (BLAST output), but found {df.shape[1]} in {path}"
        )
    df.columns = columns

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
    lineage_df["full_lineage"] = lineage_df["full_lineage"].fillna("").astype(str)
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

def filter_and_format(df, sample_name, filter_file, headers):
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

    df["ncontigs_per_sacc"] = df.groupby(["species_updated", "sacc"])["sacc"].transform("count")
    df["ncontigs_per_spp"] = df.groupby("species_updated").species_updated.transform("size")
          
    #apply scores at species level
    df["pident"] = df["pident"].astype(float)
    df["bitscore"] = df["bitscore"].astype(int)
    df["evalue"] = df["evalue"].astype(float)
    #modify how we access the cov value from the header
    headers_df = pd.read_csv(headers, sep="\t", header=None, names=["spades_headers"], dtype=str)
    headers_df["assembly_kmer_cov"] = headers_df["spades_headers"].str.extract(r'_cov_([0-9.]+)').astype(float)
    headers_df["qseqid"] = headers_df["spades_headers"].str.extract(r'^(CONTIG_[0-9]+)').astype(str)
    df = df.merge(headers_df[["qseqid", "assembly_kmer_cov"]], on="qseqid", how="left")
    df["qlen"] = df["qlen"].astype(int)
    df["alignment_length"] = df["alignment_length"].astype(int)
    df = apply_group_score(df, "species", "pident", max_pid, "pident_score_spp")
    df = apply_group_score(df, "species", "bitscore", max_bitscore, "bitscore_score_spp")
    df = apply_group_score(df, "species", "evalue", min_evalue, "evalue_score_spp")
    df = apply_group_score(df, "species", "assembly_kmer_cov", max_assembly_kmer_cov, "assembly_kmer_cov_score_spp")
    df["qcovs_score"] = global_qcovs(df)
    df = apply_group_score(df, "species", "qcovs", max_qcovs, "best_qcovs_score_spp")
    df["completeness_score"] = df["stitle"].apply(completeness_score)
    df = apply_group_score(df, "species", "qlen", max_length, "query_length_score_spp")
    df = apply_group_score(df, "species", "alignment_length", max_length, "alignment_length_score_spp")

    df["total_score_spp"] = (
            df["pident_score_spp"]
            + df["bitscore_score_spp"]
            + df["evalue_score_spp"]
            + df["assembly_kmer_cov_score_spp"]
            + df["qcovs_score"]
            + df["best_qcovs_score_spp"]
            + df["completeness_score"]
            + df["alignment_length_score_spp"]
            + df["query_length_score_spp"]
    )
    best_idx_per_spp = df.groupby("species")["total_score_spp"].idxmax()
    df["best_contig_per_sp_filter"] = df.index.isin(best_idx_per_spp)

    #apply scores at species-rna level
    df = apply_group_score(df, "species_updated", "pident", max_pid, "pident_score_spp_rna")
    df = apply_group_score(df, "species_updated", "bitscore", max_bitscore, "bitscore_score_spp_rna")
    df = apply_group_score(df, "species_updated", "evalue", min_evalue, "evalue_score_spp_rna")
    df = apply_group_score(df, "species_updated", "assembly_kmer_cov", max_assembly_kmer_cov, "assembly_kmer_cov_score_spp_rna")
    df["qcovs_score"] = global_qcovs(df)
    df = apply_group_score(df, "species_updated", "qcovs", max_qcovs, "best_qcovs_score_spp_rna")
    df["completeness_score"] = df["stitle"].apply(completeness_score)
    df = apply_group_score(df, "species_updated", "qlen", max_length, "query_length_score_spp_rna")
    df = apply_group_score(df, "species_updated", "alignment_length", max_length, "alignment_length_score_spp_rna")

    df["total_score_spp_rna"] = (
            df["pident_score_spp_rna"]
            + df["bitscore_score_spp_rna"]
            + df["evalue_score_spp_rna"]
            + df["assembly_kmer_cov_score_spp_rna"]
            + df["qcovs_score"]
            + df["best_qcovs_score_spp_rna"]
            + df["completeness_score"]
            + df["alignment_length_score_spp_rna"]
            + df["query_length_score_spp_rna"]
    )
    best_idx_per_spp_rna = df.groupby("species_updated")["total_score_spp_rna"].idxmax()
    df["best_contig_per_sp_rna_filter"] = df.index.isin(best_idx_per_spp_rna)


    #apply the same scores per accession number
    df = apply_group_score(df, "sacc", "pident", max_pid, "pident_score_sacc")
    df = apply_group_score(df, "sacc", "bitscore", max_bitscore, "bitscore_score_sacc")
    df = apply_group_score(df, "sacc", "evalue", min_evalue, "evalue_score_sacc")
    df = apply_group_score(df, "sacc", "assembly_kmer_cov", max_assembly_kmer_cov, "assembly_kmer_cov_score_sacc")
    df = apply_group_score(df, "sacc", "qcovs", max_qcovs, "best_qcovs_score_sacc")
    df = apply_group_score(df, "sacc", "qlen", max_length, "query_length_score_sacc")
    df = apply_group_score(df, "sacc", "alignment_length", max_length, "alignment_length_score_sacc")

    df["total_score_sacc"] = (
            df["pident_score_sacc"]
            + df["bitscore_score_sacc"]
            + df["evalue_score_sacc"]
            + df["assembly_kmer_cov_score_sacc"]
            + df["qcovs_score"]
            + df["best_qcovs_score_sacc"]
            + df["completeness_score"]
            + df["alignment_length_score_sacc"]
            + df["query_length_score_sacc"]
    )

    # Find the index of the best hit per species
    best_idx_per_acc = df.groupby("sacc")["total_score_sacc"].idxmax()

    #Create a new column indicating whether the row is the best hit
    df["best_contig_per_acc_filter"] = df.index.isin(best_idx_per_acc)

    # Read exclusion patterns from a file
    with open(filter_file, "r") as f:
        exclude_patterns = [line.strip() for line in f if line.strip()]

    pattern = "|".join(exclude_patterns)

    # Create a new column that indicates whether the row is filtered
    df["term_filter"] = ~(
        df["species_updated"].str.contains(pattern, case=False, na=False)
        | df["stitle"].str.contains(pattern, case=False, na=False)
        )
    # Filter out rows where cov < 5
    # Filter out rows where qcovs < 30
    #top_hits_df = top_hits_df[top_hits_df["assembly_kmer_cov"] >= 5].copy()
    #top_hits_df = top_hits_df[top_hits_df["qcovs"] >= 30].copy()    
    
    df["cov_filter"] = (
        (df["assembly_kmer_cov"] >= 1) &
        (df["qcovs"] >= 30)
        )
    
    final_columns_filt = ["sample_name", "qseqid", "sacc", "alignment_length", "evalue", "bitscore", "pident", "mismatch",
                     "gapopen", "qstart", "qend", "qlen", "sstart", "send", "slen", "sstrand",
                     "qcovhsp", "staxids", "qseq", "sseq", "qcovs", "species",
                     "species_updated", "RNA_type", "stitle", "full_lineage", "ncontigs_per_sacc", "ncontigs_per_spp", "assembly_kmer_cov", "total_score_spp", "total_score_spp_rna", "total_score_sacc", "term_filter", "cov_filter", "best_contig_per_sp_filter", "best_contig_per_sp_rna_filter", "best_contig_per_acc_filter"]
    
    return df[final_columns_filt]

def apply_group_score(df, group_col, target_col, score_func, new_col):
        """
        Apply a grouped scoring function (max, min, custom) to a column.

        Parameters
        ----------
        df : DataFrame
        group_col : str
            Column to group by.
        target_col : str
            Column to compute the score on.
        score_func : callable
            Function used inside transform(), e.g., max, min.
        new_col : str
            Output column name.
        """
        df[new_col] = df.groupby(group_col)[target_col].transform(score_func)
        return df
  

def max_naccs(series):
    max_val =  series.max()
    return (series == max_val).astype(int)

def max_pid(series):
    max_val = series.max()
    return (series == max_val).astype(int)

def max_bitscore(series):
    max_val = series.max()
    return (series == max_val).astype(int) * 2

def max_length(series):
    max_val = series.max()
    return (series == max_val).astype(int) * 2

def min_evalue(series):
    min_val = series.min()
    return (series == min_val).astype(int) * 2

def global_qcovs(df):
    return (df["qcovs"] > 75).astype(int)

def max_qcovs(series):
    max_val = series.max()
    return (series == max_val).astype(int) * 2

def max_assembly_kmer_cov(series):
    max_val = series.max()
    return (series == max_val).astype(int) * 2

def completeness_score(x):
    if "complete sequence" in str(x):
        return 3
    elif "complete genome" in str(x):
        return 3
    elif "polyprotein gene complete cds" in str(x):
        return 3
    elif "polyprotein 1 gene complete cds" in str(x):
        return 3
    elif "polyprotein 2 gene complete cds" in str(x):
        return 3
    elif "nearly complete sequence" in str(x):
        return -3
    if "partial" in str(x):
        return -3
    elif "polymerase protein" in str(x):
        return -3
    elif "RNA-dependent RNA polymerase" in str(x):
        return -3
    else:   
        return 0
    
def main():
    args = parse_arguments()
    blastn_results_path = args.blastn_results
    sample_name = args.sample_name
    tk_db_dir = args.taxonkit_database_dir
    filter_file = args.filter
    headers = args.assembly_headers

    if not os.path.isfile(blastn_results_path):
        raise FileNotFoundError(f"{blastn_results_path} does not exist.")

    blastn_results = load_blast_results(blastn_results_path)
    enriched_dfs = enrich_with_taxonomy(blastn_results, tk_db_dir)
    merged_df = merge_taxonomy(enriched_dfs)
    final_df = filter_and_format(merged_df, sample_name, filter_file, headers)

    out_file = os.path.basename(args.blastn_results).replace("_blastn.txt", "_megablast_top_viral_hits.txt")
    final_df.to_csv(out_file, sep="\t", index=False)
    print(f"Virus-only results saved to {out_file}")

if __name__ == "__main__":
    main()

#singularity exec  -B /work/daff_viral_rnaseq/rnaspades/BR_Guar_18 docker://quay.io/biocontainers/pytaxonkit:0.9.1--pyhdfd78af_1 /work/daff_viral_rnaseq/rnaspades/BR_Guar_18/filter_blast.py --blastn_results BR_Guar_18_blastn_header.txt --sample_name BR_Guar_18 --mode ncbi --taxonkit_database_dir ~/.taxonkit
