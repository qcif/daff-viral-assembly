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
    return parser.parse_args()

def load_blast_results(path):
    """Load BLASTn results based on mode."""
    if not os.path.isfile(path):
        raise FileNotFoundError(f"BLASTn results file not found: {path}")
    # Define column names as a string
    header = ("qseqid\tsgi\tsacc\tlength\tpident\tmismatch\tgapopen\tqstart\tqend\tqlen\t"
              "sstart\tsend\tslen\tsstrand\tevalue\tbitscore\tqcovhsp\tstitle\tstaxids\tqseq\t"
              "sseq\tsseqid\tqcovs\tqframe\tsframe\n")
    
    # Create a temporary file with header + original content
    with tempfile.NamedTemporaryFile(mode='w+', delete=False) as tmp_file:
        tmp_file.write(header)
        with open(path, 'r') as original_file:
            shutil.copyfileobj(original_file, tmp_file)
        tmp_path = tmp_file.name


    # Define columns and dtypes
    columns = header.strip().split('\t')
    #columns = ["qseqid", "sgi", "sacc", "length", "nident", "pident", "mismatch", "gaps", "gapopen", "qstart",
    #           "qend", "qlen", "sstart", "send", "slen", "sstrand", "evalue", "bitscore", "qcovhsp", "stitle",
    #           "staxids", "qseq", "sseq", "sseqid", "qcovs", "qframe", "sframe"]

    dtype = {
        "qseqid": 'str', "sgi": 'str', "sacc": 'str', "length": 'int64', "pident": 'float64', "mismatch": 'int64',
        "gapopen": 'int64', "qstart": 'int64', "qend": 'int64', "qlen": 'int64', "sstart": 'int64', 
        "send": 'int64', "slen": 'int64', "sstrand": 'str', "evalue": 'float64', "bitscore": 'float64', 
        "qcovhsp": 'int64', "stitle": 'str', "staxids": 'str', "qseq": 'str', "sseq": 'str', "sseqid": 'str', 
        "qcovs": 'int64', "qframe": 'int64', "sframe": 'int64'
    }
    # Load DataFrame
    df = pd.read_csv(tmp_path, sep='\t', usecols=columns, dtype=dtype)
    #df = pd.read_csv(path, sep="\t", header=0, usecols=columns, dtype=dtype)
    df["staxids"] = pd.to_numeric(df["staxids"].str.split(";").str[0], errors='coerce').fillna(0).astype(int)
    #top_hit = df.drop_duplicates(subset=["qseqid"], keep="first").copy()

    return df

def enrich_with_taxonomy(df, taxonkit_dir):
    """Add taxonomy information to the DataFrame."""
    staxids_l = df["staxids"].unique().tolist()

    lineage_df = pytaxonkit.lineage(staxids_l, data_dir=taxonkit_dir)[['TaxID', 'FullLineage']]
    lineage_df.columns = ["staxids", "FullLineage"]
    lineage_df["staxids"] = lineage_df["staxids"].astype(int)
    lineage_df["FullLineage"] = lineage_df["FullLineage"].str.lower().str.replace(" ", "_", regex=False)
    lineage_df["broad_taxonomic_category"] = np.where(
        lineage_df["FullLineage"].str.contains("viruses;"),
        "virus",
        np.where(
            lineage_df["FullLineage"].str.contains("viroids;"),
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

def filter_and_format(df, sample_name, filter):
    """Filter only viral hits and format."""
    df.insert(0, "sample_name", sample_name)
    df = df[~df["species"].str.contains("synthetic construct", na=False)]
    df = df[~df["species"].str.contains("Expression vector", na=False)]
    #df = df.drop_duplicates(subset=["qseqid"], keep="first").copy()
    # Then drop duplicates — keep first valid hit per qseqid
    df = df.sort_values(by=["qseqid", "bitscore"], ascending=[True, False])  # Optional: Sort by best score
    df = df.drop_duplicates(subset=["qseqid"], keep="first").copy()
    

    # Filter only viral entries
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

    #df["Species_updated"] = df[["species", "RNA_type"]].agg(" ".join, axis=1)
    df["Species_updated"] = df[["species", "RNA_type"]].fillna("").astype(str).agg(" ".join, axis=1)
    df["Species_updated"] = df["Species_updated"].str.replace("RNA1 RNA1", "RNA1", regex=False)
    df["Species_updated"] = df["Species_updated"].str.replace("RNA2 RNA2", "RNA2", regex=False)
    df["Species_updated"] = df["Species_updated"].str.replace("RNA3 RNA3", "RNA3", regex=False)
    df["Species_updated"] = df["Species_updated"].str.rstrip()

    df["naccs"] = df.groupby(["Species_updated", "sacc"])["sacc"].transform("count")
    df["naccs_score"] = df.groupby("Species_updated")["naccs"].transform(max_naccs)

    df["pident"] = df["pident"].astype(float)
    df["pident_score"] = df.groupby("Species_updated")["pident"].transform(max_pid)

    df["bitscore"] = df["bitscore"].astype(int)
    df["bitscore_score"] = df.groupby("Species_updated")["bitscore"].transform(max_bitscore)

    df["length"] = df["length"].astype(int)
    df["length_score"] = df.groupby("Species_updated")["length"].transform(max_length)

    df["cov"] = df["qseqid"].str.extract(r'_cov_([0-9.]+)').astype(float)
    df["cov_score"] = df.groupby("Species_updated")["cov"].transform(max_cov)

    df["evalue"] = df["evalue"].astype(float)
    df["evalue_score"] = df.groupby("Species_updated")["evalue"].transform(min_evalue)
    # assign a score of 1 if qcovs > 75, else 0
    df["global_qcovs_score"] = global_qcovs(df)
    df["qcovs_score"] = df.groupby("Species_updated")["qcovs"].transform(max_qcovs)
    df["completeness_score"] = df["stitle"].apply(completeness_score)
    df["total_score"] = df["naccs_score"] + df["pident_score"] + df["bitscore_score"] + df["length_score"] + df["evalue_score"] + df["global_qcovs_score"] + df["qcovs_score"] + df["completeness_score"] + df["cov_score"].astype(int)

    
    final_columns = ["sample_name", "qseqid", "sgi", "sacc", "length", "pident", "mismatch",
                     "gapopen", "qstart", "qend", "qlen", "sstart", "send", "slen", "sstrand", "evalue", "bitscore",
                     "qcovhsp", "stitle", "staxids", "qseq", "sseq", "sseqid", "qcovs", "qframe", "sframe",
                     "species", "Species_updated", "broad_taxonomic_category", "FullLineage", "naccs", "naccs_score", "pident_score",
                     "bitscore_score", "evalue_score", "global_qcovs_score", "cov_score", "qcovs_score", "completeness_score", "length_score", "total_score"]
    

    top_hits_df = df.loc[df.groupby(["Species_updated"])["total_score"].idxmax()].copy()
    # Read exclusion patterns from a file
    with open(filter, "r") as f:
        exclude_patterns = [line.strip() for line in f if line.strip()]

    pattern = "|".join(exclude_patterns)

    top_hits_df = top_hits_df[~top_hits_df["Species_updated"].str.contains(pattern, case=False, na=False)]
    top_hits_df = top_hits_df[~top_hits_df["stitle"].str.contains(pattern, case=False, na=False)]
    # Extract cov value using regex and convert to float
    top_hits_df["cov"] = top_hits_df["qseqid"].str.extract(r'_cov_([0-9.]+)').astype(float)

    # Filter out rows where cov < 5
    top_hits_df = top_hits_df[top_hits_df["cov"] >= 1].copy()
    top_hits_df = top_hits_df[top_hits_df["qcovs"] >= 30].copy()    
    #
    #exclude_patterns = ["phage", "tick virus", "Sclerotinia", "Plasmopara", "Botrytis cinerea", "Erysiphales","Erysiphe"]
    #pattern = "|".join(exclude_patterns)  # create regex pattern: 'phage|tick virus|Sclerotinia|Plasmopara'
    #top_hits_df = top_hits_df[~top_hits_df["species"].str.contains(pattern, case=False, na=False)]

    return df[final_columns], top_hits_df

def max_naccs(series):
    max_val = max_val = series.max()
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
    min_val = series.min()
    return (series == min_val).astype(int) * 2

def max_cov(series):
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
    filter = args.filter

    if not os.path.isfile(blastn_results_path):
        raise FileNotFoundError(f"{blastn_results_path} does not exist.")

    blastn_results = load_blast_results(blastn_results_path)
    enriched_dfs = enrich_with_taxonomy(blastn_results, tk_db_dir)
    merged_df = merge_taxonomy(enriched_dfs)
    final_df, filtered_df = filter_and_format(merged_df, sample_name, filter)

    out_file = os.path.basename(args.blastn_results).replace("_blastn.txt", "_megablast_top_viral_hits.txt")
    final_df.to_csv(out_file, sep="\t", index=False)
    print(f"Virus-only results saved to {out_file}")
    
    out_file2 = os.path.basename(args.blastn_results).replace("_blastn.txt", "_megablast_top_viral_hits_filtered.txt")
    filtered_df.to_csv(out_file2, sep="\t", index=False)
    print(f"Virus-only filtered results saved to {out_file2}")

if __name__ == "__main__":
    main()

#singularity exec  -B /work/daff_viral_rnaseq/rnaspades/BR_Guar_18 docker://quay.io/biocontainers/pytaxonkit:0.9.1--pyhdfd78af_1 /work/daff_viral_rnaseq/rnaspades/BR_Guar_18/filter_blast.py --blastn_results BR_Guar_18_blastn_header.txt --sample_name BR_Guar_18 --mode ncbi --taxonkit_database_dir ~/.taxonkit
