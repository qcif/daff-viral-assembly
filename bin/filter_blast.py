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
    #header = ("qseqid\tsgi\tsacc\talignment_length\tpident\tmismatch\tgapopen\tqstart\tqend\tqlen\t"
    #          "sstart\tsend\tslen\tsstrand\tevalue\tbitscore\tqcovhsp\tstitle\tstaxids\tqseq\t"
    #          "sseq\tsseqid\tqcovs\tqframe\tsframe\n")
    # Define expected header (based on your BLAST output fields)
    columns = [
        "qseqid", "sgi", "sacc", "alignment_length", "pident", "mismatch", "gapopen",
        "qstart", "qend", "qlen", "sstart", "send", "slen", "sstrand", "evalue",
        "bitscore", "qcovhsp", "stitle", "staxids", "qseq", "sseq", "sseqid",
        "qcovs", "qframe", "sframe"
    ]
    # Create a temporary file with header + original content
    #with tempfile.NamedTemporaryFile(mode='w+', delete=False) as tmp_file:
    #    tmp_file.write(header)
    #    with open(path, 'r') as original_file:
    #        shutil.copyfileobj(original_file, tmp_file)
    #    tmp_path = tmp_file.name
    df = pd.read_csv(path, sep="\t", header=None, dtype=str)
    if df.shape[1] != len(columns):
        raise ValueError(
            f"Expected {len(columns)} columns (BLAST output), but found {df.shape[1]} in {path}"
        )
    df.columns = columns

    # Define columns and dtypes
    #columns = header.strip().split('\t')
    #columns = ["qseqid", "sgi", "sacc", "length", "nident", "pident", "mismatch", "gaps", "gapopen", "qstart",
    #           "qend", "qlen", "sstart", "send", "slen", "sstrand", "evalue", "bitscore", "qcovhsp", "stitle",
    #           "staxids", "qseq", "sseq", "sseqid", "qcovs", "qframe", "sframe"]

    dtype = {
        "qseqid": 'str', "sgi": 'str', "sacc": 'str', "alignment_length": 'int64', "pident": 'float64', "mismatch": 'int64',
        "gapopen": 'int64', "qstart": 'int64', "qend": 'int64', "qlen": 'int64', "sstart": 'int64', 
        "send": 'int64', "slen": 'int64', "sstrand": 'str', "evalue": 'float64', "bitscore": 'float64', 
        "qcovhsp": 'int64', "stitle": 'str', "staxids": 'str', "qseq": 'str', "sseq": 'str', "sseqid": 'str', 
        "qcovs": 'int64', "qframe": 'int64', "sframe": 'int64'
    }
    # Load DataFrame
    #df = pd.read_csv(tmp_path, sep='\t', usecols=columns, dtype=dtype)
    #df = pd.read_csv(path, sep="\t", header=0, usecols=columns, dtype=dtype)
    for col, dtype_ in dtype.items():
        if col in df.columns:
            if dtype_ in ("int64", "float64"):
                df[col] = pd.to_numeric(df[col], errors="coerce").astype(dtype_)
            else:
                df[col] = df[col].astype(str)

    df["staxids"] = pd.to_numeric(df["staxids"].str.split(";").str[0], errors='coerce').fillna(0).astype(int)
    #top_hit = df.drop_duplicates(subset=["qseqid"], keep="first").copy()

    return df

def enrich_with_taxonomy(df, taxonkit_dir):
    """Add taxonomy information to the DataFrame."""
    staxids_l = df["staxids"].unique().tolist()
    #lineage_df = pytaxonkit.lineage(staxids_l, data_dir=taxonkit_dir)
    #lineage_df.rename(columns=str.lower, inplace=True)      
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

def filter_and_format(df, sample_name, filter_file):
    """Filter only viral hits and format."""
    df.insert(0, "sample_name", sample_name)
    df = df[~df["species"].str.contains("synthetic construct", na=False)]
    df = df[~df["species"].str.contains("Expression vector", na=False)]
    #df = df.drop_duplicates(subset=["qseqid"], keep="first").copy()
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

    #df["species_updated"] = df[["species", "RNA_type"]].agg(" ".join, axis=1)
    df["species_updated"] = df[["species", "RNA_type"]].fillna("").astype(str).agg(" ".join, axis=1)
    df["species_updated"] = df["species_updated"].str.replace("RNA1 RNA1", "RNA1", regex=False)
    df["species_updated"] = df["species_updated"].str.replace("RNA2 RNA2", "RNA2", regex=False)
    df["species_updated"] = df["species_updated"].str.replace("RNA3 RNA3", "RNA3", regex=False)
    df["species_updated"] = df["species_updated"].str.rstrip()

    df["ncontigs"] = df.groupby(["species_updated", "sacc"])["sacc"].transform("count")
    df["ncontigs_score"] = df.groupby("species_updated")["ncontigs"].transform(max_naccs)

    df["pident"] = df["pident"].astype(float)
    #df["pident"] = pd.to_numeric(df["pident"], errors="coerce")
    df["pident_score"] = df.groupby("species_updated")["pident"].transform(max_pid)

    df["bitscore"] = df["bitscore"].astype(int)
    #df["bitscore"] = pd.to_numeric(df["bitscore"], errors="coerce")
    df["bitscore_score"] = df.groupby("species_updated")["bitscore"].transform(max_bitscore)

    df["evalue"] = df["evalue"].astype(float)
    #df["evalue"] = pd.to_numeric(df["evalue"], errors="coerce")
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
    #df["alignment_length"] = pd.to_numeric(df["alignment_length"], errors="coerce")
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

    #top_hits_df = df.loc[df.groupby(["species_updated"])["total_score"].idxmax()].copy()


    # Find the index of the best hit per species
    best_idx = df.groupby("species_updated")["total_score"].idxmax()

    #Create a new column indicating whether the row is the best hit
    df["best_contig_per_sp_filter"] = df.index.isin(best_idx)


    # Read exclusion patterns from a file
    with open(filter_file, "r") as f:
        exclude_patterns = [line.strip() for line in f if line.strip()]

    pattern = "|".join(exclude_patterns)

    #top_hits_df = top_hits_df[~top_hits_df["species_updated"].str.contains(pattern, case=False, na=False)]
    #top_hits_df = top_hits_df[~top_hits_df["stitle"].str.contains(pattern, case=False, na=False)]
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
        (df["assembly_kmer_cov"] >= 5) &
        (df["qcovs"] >= 30)
        )

    #exclude_patterns = ["phage", "tick virus", "Sclerotinia", "Plasmopara", "Botrytis cinerea", "Erysiphales","Erysiphe"]
    #pattern = "|".join(exclude_patterns)  # create regex pattern: 'phage|tick virus|Sclerotinia|Plasmopara'
    #top_hits_df = top_hits_df[~top_hits_df["species"].str.contains(pattern, case=False, na=False)]

    #final_columns = ["sample_name", "qseqid", "sacc", "alignment_length", "evalue", "bitscore", "pident", "mismatch",
    #                 "gapopen", "qstart", "qend", "qlen", "sstart", "send", "slen", "sstrand", 
    #                 "qcovhsp", "staxids", "qseq", "sseq", "qcovs", 
    #                 "species_updated", "RNA_type", "stitle", "full_lineage", "ncontigs", 
    #                 "ncontigs_score", "pident_score", "bitscore_score", "evalue_score", "assembly_kmer_cov_score", 
    #                 "qcovs_score", "best_qcovs_score", "completeness_score", "alignment_length_score", "query_length_score", "total_score"]
    
    final_columns_filt = ["sample_name", "qseqid", "sacc", "alignment_length", "evalue", "bitscore", "pident", "mismatch",
                     "gapopen", "qstart", "qend", "qlen", "sstart", "send", "slen", "sstrand",
                     "qcovhsp", "staxids", "qseq", "sseq", "qcovs", 
                     "species_updated", "RNA_type", "stitle", "full_lineage", "ncontigs", "total_score", "term_filter", "cov_filter", "best_contig_per_sp_filter"]
    
    #columns before filtering and re-ordering:
    # sample_name	qseqid	sgi	sacc	alignment_length	pident	mismatch	gapopen	qstart	qend	qlen	sstart	send	slen	
    # sstrand	evalue	bitscore	qcovhsp	stitle	staxids	qseq	sseq	sseqid	qcovs	qframe	sframe	species	full_lineage	
    # broad_taxonomic_category	RNA_type	species_updated	ncontigs	ncontigs_score	pident_score	bitscore_score	alignment_length_score	
    # cov	cov_score	evalue_score	global_qcovs_score	qcovs_score	completeness_score	total_score

    #return df[final_columns], top_hits_df[final_columns_filt]
    return df[final_columns_filt]
    

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

    if not os.path.isfile(blastn_results_path):
        raise FileNotFoundError(f"{blastn_results_path} does not exist.")

    blastn_results = load_blast_results(blastn_results_path)
    enriched_dfs = enrich_with_taxonomy(blastn_results, tk_db_dir)
    merged_df = merge_taxonomy(enriched_dfs)
    #final_df, filtered_df = filter_and_format(merged_df, sample_name, filter_file)
    final_df = filter_and_format(merged_df, sample_name, filter_file)

    out_file = os.path.basename(args.blastn_results).replace("_blastn.txt", "_megablast_top_viral_hits.txt")
    final_df.to_csv(out_file, sep="\t", index=False)
    print(f"Virus-only results saved to {out_file}")
    
    #out_file2 = os.path.basename(args.blastn_results).replace("_blastn.txt", "_megablast_top_viral_hits_filtered.txt")
    #filtered_df.to_csv(out_file2, sep="\t", index=False)
    #print(f"Virus-only filtered results saved to {out_file2}")

if __name__ == "__main__":
    main()

#singularity exec  -B /work/daff_viral_rnaseq/rnaspades/BR_Guar_18 docker://quay.io/biocontainers/pytaxonkit:0.9.1--pyhdfd78af_1 /work/daff_viral_rnaseq/rnaspades/BR_Guar_18/filter_blast.py --blastn_results BR_Guar_18_blastn_header.txt --sample_name BR_Guar_18 --mode ncbi --taxonkit_database_dir ~/.taxonkit
