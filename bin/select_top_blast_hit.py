#!/usr/bin/env python

import pandas as pd
import argparse
import os
from functools import reduce
import pytaxonkit
import numpy as np
import re

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Load and enrich BLASTn results.")
    parser.add_argument("--blastn_results", required=True, type=str)
    parser.add_argument("--sample_name", required=True, type=str)
    parser.add_argument("--target_organism", required=True, type=str)
    parser.add_argument("--taxonkit_database_dir", required=True, type=str)
    return parser.parse_args()

def load_blast_results(path):
    """Load BLASTn results based on mode."""
    if not os.path.isfile(path):
        raise FileNotFoundError(f"BLASTn results file not found: {path}")

    columns = ["qseqid", "sgi", "sacc", "length", "nident", "pident", "mismatch", "gaps", "gapopen", "qstart",
                "qend", "qlen", "sstart", "send", "slen", "sstrand", "evalue", "bitscore", "qcovhsp", "stitle",
                "staxids", "qseq", "sseq", "sseqid", "qcovs", "qframe", "sframe"]

    dtype = {
        "qseqid": 'str', "sgi": 'str', "sacc": 'str', "length": 'int64', "nident": 'int64',
        "pident": 'float64', "mismatch": 'int64', "gaps": 'int64', "gapopen": 'int64', "qstart": 'int64',
        "qend": 'int64', "qlen": 'int64', "sstart": 'int64', "send": 'int64', "slen": 'int64', "sstrand": 'str',
        "evalue": 'float64', "bitscore": 'float64', "qcovhsp": 'int64', "stitle": 'str', "staxids": 'str',
        "qseq": 'str', "sseq": 'str', "sseqid": 'str', "qcovs": 'int64', "qframe": 'int64', "sframe": 'int64'
    }

    df = pd.read_csv(path, sep="\t", header=0, usecols=columns, dtype=dtype)
    df["staxids"] = pd.to_numeric(df["staxids"].str.split(";").str[0], errors='coerce').fillna(0).astype(int)
    top_hit = df.drop_duplicates(subset=["qseqid"], keep="first").copy()

    return top_hit

def enrich_with_taxonomy(df, taxonkit_dir):
    """Add taxonomy information to the dataFrame."""
    #retain unique staxids
    staxids_l = df["staxids"].unique().tolist()

    lineage_df = pytaxonkit.lineage(staxids_l, data_dir=taxonkit_dir)[['TaxID', 'FullLineage']]
    lineage_df.columns = ["staxids", "FullLineage"]
    lineage_df["staxids"] = lineage_df["staxids"].astype(int)
    lineage_df["FullLineage"] = lineage_df["FullLineage"].str.lower().str.replace(" ", "_", regex=False)
    lineage_df["broad_taxonomic_category"] = np.where(
        lineage_df["FullLineage"].str.contains("viruses;"),
        "virus",
        np.where(
            lineage_df["FullLineage"].str.contains(";candidatus_phytoplasma;"),
            "bacteria;phytoplasma",
            np.where(
            (lineage_df["FullLineage"].str.contains(";bacteria;")) &
            (~lineage_df["FullLineage"].str.contains(";candidatus_phytoplasma;")),
            "bacteria;other",
                np.where(
                    lineage_df["FullLineage"].str.contains(";archaea;"),
                    "archaea",
                    np.where(
                        lineage_df["FullLineage"].str.contains(";erysiphaceae;"),
                        "eukaryota;fungi;powdery_mildew",
                        np.where(
                            (lineage_df["FullLineage"].str.contains(";fungi;")) &
                            (~lineage_df["FullLineage"].str.contains(";erysiphaceae;")),
                            "eukaryota;fungi;other",
                            np.where(
                                (lineage_df["FullLineage"].str.contains(";deuterostomia;")),
                                "eukaryota;deuterostomia",
                                np.where(
                                    (lineage_df["FullLineage"].str.contains(";protostomia;")),
                                    "eukaryota;protostomia",
                                    np.where(
                                        (lineage_df["FullLineage"].str.contains(";eukaryota;")) &
                                        (~lineage_df["FullLineage"].str.contains(";fungi;")) &
                                        (~lineage_df["FullLineage"].str.contains(";deuterostomia;")) &
                                        (~lineage_df["FullLineage"].str.contains(";protostomia;")),
                                        "eukaryota;other",
                                        "other"
                                    )
                                )
                            )
                        )
                    )
                )
            )
        )
    )
    print(lineage_df)

    names_df = pytaxonkit.name(staxids_l, data_dir=taxonkit_dir)[['TaxID', 'Name']]
    names_df.columns = ["staxids", "species"]
    names_df["staxids"] = names_df["staxids"].astype(int)

    return [df, names_df, lineage_df]

def merge_taxonomy(dfs):
    """Merge taxonomy-enriched data."""
    return reduce(lambda left, right: pd.merge(left, right, on="staxids", how="outer"), dfs)

def filter_and_format(df, sample_name, target_organism):
    """Final formatting, filtering, and matching."""
    df.insert(0, "sample_name", sample_name)
    df = df[~df["species"].str.contains("synthetic construct", na=False)]

    # Enable to have a list of target organisms.
    # Normalise target_organism to a list of lowercase, underscored strings
    if '|' in target_organism:
        target_organism = [s.strip().strip("'\"") for s in target_organism.split('|')]
        target_organisms = [org.lower().replace(" ", "_") for org in target_organism]
    else:
        target_organisms = [target_organism.lower().replace(" ", "_")]

    # Create a combined mask for all target organisms
    mask = False
    # Check if the target organism is in the broad taxonomic category or in the full lineage
    for org in target_organisms:
        pattern = rf"\b{re.escape(org)}\b"  # \b ensures word boundary match

        exact_broad_match = df["broad_taxonomic_category"].str.lower().str.contains(pattern, na=False, regex=True)
        exact_lineage_match = df["FullLineage"].apply(lambda x: org in [s.strip().lower() for s in str(x).split(';')]
        )
        mask |= exact_broad_match | exact_lineage_match
    df["target_organism_match"] = np.where(mask, "Y", "N")

    
    final_columns = ["sample_name", "qseqid", "sgi", "sacc", "length", "nident", "pident", "mismatch", "gaps",
                     "gapopen", "qstart", "qend", "qlen", "sstart", "send", "slen", "sstrand", "evalue", "bitscore",
                     "qcovhsp", "stitle", "staxids", "qseq", "sseq", "sseqid", "qcovs", "qframe", "sframe",
                     "species", "broad_taxonomic_category", "FullLineage", "target_organism_match"]

    return df[final_columns]

def main():
    args = parse_arguments()
    blastn_results_path = args.blastn_results
    sample_name = args.sample_name
    target_organism = args.target_organism
    tk_db_dir = args.taxonkit_database_dir

    if not os.path.isfile(blastn_results_path):
        raise FileNotFoundError(f"{blastn_results_path} does not exist.")
    #elif len(blastn_results_path) == 0:
    #    print("DataFrame is empty!")
    #    out_file = open(os.path.basename(args.blastn_results).replace("_top_10_hits.txt", "_top_hits_tmp.txt"))
    #    out_file.write("sample_name\tqseqid\tsgi\tsacc\tlength\tnident\tpident\tmismatch\tgaps\tgapopen\tqstart\tqend\tqlen\tsstart\tsend\tslen\tsstrand\tevalue\tbitscore\tqcovhsp\tstitle\tstaxids\tqseq\tsseq\tsseqid\tqcovs\tqframe\tsframe\tspecies\tbroad_taxonomic_category\tFullLineage\ttarget_organism_match")
    #    out_file.close()
    #    exit ()

    blastn_results = load_blast_results(blastn_results_path)
    enriched_dfs = enrich_with_taxonomy(blastn_results, tk_db_dir)
    merged_df = merge_taxonomy(enriched_dfs)
    final_df = filter_and_format(merged_df, sample_name, target_organism)
    out_file = os.path.basename(args.blastn_results).replace("_top_10_hits.txt", "_top_hits_tmp.txt")
    final_df.to_csv(out_file, sep="\t", index=False)
    print(f"Results saved to {out_file}")

if __name__ == "__main__":
    main()