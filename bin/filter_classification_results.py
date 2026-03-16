#!/usr/bin/env python
import pandas as pd
import argparse
import os
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
    parser.add_argument("--outdir", default=".", help="output directory")
    return parser.parse_args()

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

#def categorise_kaiju(row):
#    #lineage = str(row.get("full_lineage", "")).lower()
#    viral_category = classify_viral_resolution(row)
#    if viral_category:
#        return viral_category
#    name = str(row.get("taxon_name", "")).lower()
##    if "unclassified" in name:
#        return "unclassified"
#    return "other non-viral"
def categorise_kaiju(row):
    resolution = classify_taxonomic_resolution(row)

    if resolution:
        return resolution

    name = str(row.get("taxon_name", "")).lower()

    if "unclassified" in name:
        return "unclassified"

    return "other"


#def categorise_bracken(row):
#    viral_category = classify_viral_resolution(row)
#    if viral_category:
#        return viral_category
#    name = str(row.get("taxon_name", "")).lower()
#    if "unclassified" in name:
##        return "unclassified"
#    return "other non-viral"

def categorise_bracken(row):
    resolution = classify_taxonomic_resolution(row)

    if resolution:
        return resolution

    name = str(row.get("taxon_name", "")).lower()

    if "unclassified" in name:
        return "unclassified"

    return "other"


#def classify_viral_resolution(row):
#    ranks = str(row.get("full_lineage_ranks", "")).lower()
#    lineage = str(row.get("full_lineage", "")).lower()
#    if not ranks or ("viruses" not in lineage and "viroids" not in lineage):
#        return None
#    rank_list = [r.strip() for r in ranks.split(";") if r.strip()]

    # Find deepest informative viral rank
#    if "species" in rank_list:
#        return "viral"
#    if "genus" in rank_list:
#        return "genus_unclassified"
#    if "family" in rank_list:
#        return "family_unclassified"
#    if "order" in rank_list:
#         return "order_unclassified"
#     if "class" in rank_list:
#        return "class_unclassified"
#    return "viral_high_level"
def classify_taxonomic_resolution(row):

    ranks = str(row.get("full_lineage_ranks", "")).lower()

    if not ranks:
        return "unclassified"

    rank_list = [r.strip() for r in ranks.split(";") if r.strip()]

    if "species" in rank_list:
        return "species"

    if "genus" in rank_list:
        return "genus_unclassified"

    if "family" in rank_list:
        return "family_unclassified"

    if "order" in rank_list:
        return "order_unclassified"

    if "class" in rank_list:
        return "class_unclassified"

    if "phylum" in rank_list:
        return "phylum_unclassified"

    if "kingdom" in rank_list:
        return "kingdom_unclassified"

    return "high_level_unclassified"

def classify_broad_category(row):

    lineage = str(row.get("full_lineage", "")).lower()

    if "viruses" in lineage or "viroids" in lineage:
        return "viral"

    if "unclassified" in str(row.get("taxon_name", "")).lower():
        return "unclassified_top"

    return "other"

def build_taxonomy_maps(taxids, tk_db_dir):
    lin = lineage(taxids, data_dir=tk_db_dir)
    lin_df = pd.DataFrame(lin)
    lin_df = lin_df.rename(columns={
        "TaxID": "taxid",
        "FullLineageTaxIDs": "full_lineage_taxids",
        "FullLineage": "full_lineage_names",
        "FullLineageRanks": "full_lineage_ranks",
        "Name": "name"
    })
    lin_df["full_lineage_names"] = lin_df["full_lineage_names"].fillna("NA").astype(str)
    lin_df["full_lineage_ranks"] = lin_df["full_lineage_ranks"].fillna("NA").astype(str)
    lin_df["taxid"] = pd.to_numeric(lin_df["taxid"], errors="coerce").fillna(0).astype(int)
    return (
        dict(zip(lin_df["taxid"], lin_df["full_lineage_names"])),
        dict(zip(lin_df["taxid"], lin_df["full_lineage_ranks"])),
        dict(zip(lin_df["taxid"], lin_df["full_lineage_taxids"]))
    )
def filter_viral_rows(df, pc_threshold):

    cat = df["broad_category"].astype(str).str.strip().str.lower()

    viral_labels = {
        "viral",
        "genus_unclassified",
        "family_unclassified",
        "order_unclassified",
        "class_unclassified",
        "viral_high_level"
    }

    is_viral_taxon = cat.isin(viral_labels)

    return df[
        ((is_viral_taxon) | (df["pc_reads"].astype(float) >= pc_threshold))
        & (df["reads"].astype(int) > 0)
    ].copy()

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
    #print(f"Unique taxids in Kaiju results: {len(unique_taxids)}")

    # Query lineage (just pass taxid list)
    # lineage(unique_taxids, data_dir=tk_db_dir)
    
    # Convert to DataFrame
    #lin_df = pd.DataFrame(lin)
    #print(lin_df.columns)
    #print(lin_df.iloc[0])
    # Keep relevant columns and rename
    #lin_df = lin_df.rename(columns={
    #    "TaxID": "taxid",
    #    "FullLineageTaxIDs": "full_lineage_taxids",
    #    "FullLineage": "full_lineage_names",
    #    "FullLineageRanks": "full_lineage_ranks",
    #    "Name": "name"
    #})
    taxid_to_lineage_names, taxid_to_lineage_ranks, taxid_to_lineage = \
    build_taxonomy_maps(unique_taxids, tk_db_dir)


    # Fill NaNs with 0 and convert to int
    #lin_df["taxid"] = pd.to_numeric(lin_df["taxid"], errors="coerce").fillna(0).astype(int)

    # Build mapping
    #taxid_to_lineage = {int(k): str(v) for k, v in zip(lin_df["taxid"], lin_df["full_lineage_taxids"])}

    #taxid_to_lineage_names = dict(zip(
    #    lin_df["taxid"],
    #    lin_df["full_lineage_names"]
    #))
    #taxid_to_lineage_ranks = dict(zip(
    #    lin_df["taxid"],
    #    lin_df["full_lineage_ranks"]
    #))

    # Viral check
    df["taxon_id"] = pd.to_numeric(df["taxon_id"], errors="coerce").fillna(0).astype(int)

    df["full_lineage"] = df["taxon_id"].apply(lambda x: taxid_to_lineage_names.get(int(x), ""))
    df["full_lineage_ranks"] = df["taxon_id"].apply(
        lambda x: taxid_to_lineage_ranks.get(int(x), "")
    )
    df["resolution_level"] = df.apply(classify_taxonomic_resolution, axis=1)
    df["broad_category"] = df.apply(classify_broad_category, axis=1)

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


    # Ensure reads is numeric, then filter detections with >= 10 reads
    df["reads"] = pd.to_numeric(df["reads"], errors="coerce").fillna(0).astype(int)

    df_subset = df[['taxon_name', 'taxon_id', 'full_lineage', 'full_lineage_ranks',
                    'broad_category', 'resolution_level','reads', 'pc_reads', 'term_filter', 'cov_filter']]
    
    # Apply min-reads filter ONLY to non-viral taxa
    # Keep all viral/unclassified-virus rows regardless of reads,
    # but require reads >= 10 for everything else.
    df_subset = filter_viral_rows(df_subset, 0.002)
    #cat = df_subset["broad_category"].astype(str).str.strip().str.lower()

    #viral_labels = {
    #    "viral",
    #    "genus_unclassified",
    #    "family_unclassified",
    #    "order_unclassified",
    #    "class_unclassified",
    #    "viral_high_level"
    #}
    #is_viral_taxon = cat.isin(viral_labels)

    #df_subset = df_subset[is_viral_taxon | (df_subset["pc_reads"].astype(float) >= 0.002)].copy()
    out_kaiju = os.path.join(args.outdir, f"{sample_name}_kaiju_summary.txt")
    df_subset.to_csv(out_kaiju, sep="\t", index=False)
    
    # ------------------------------------
    # Load KRAKEN/BRAKEN summary table
    # ------------------------------------
    br = pd.read_csv(bracken_path, sep="\t", dtype=str)
    br = br.rename(columns={
        "taxonomy_id": "taxon_id"})
    unique_taxids_br = br.loc[br["taxon_id"] != 0, "taxon_id"].unique().tolist()
    #lin2 = lineage(unique_taxids_br, data_dir=tk_db_dir)
    #lin2_df = pd.DataFrame(lin2)
    #lin2_df = lin2_df.rename(columns={
    #    "TaxID": "taxid",
    #    "FullLineageTaxIDs": "full_lineage_taxids",
    #    "FullLineage": "full_lineage_names",
    #    "FullLineageRanks": "full_lineage_ranks",
    #    "Name": "name"
    #})
    #lin2_df["full_lineage_ranks"] = lin2_df["full_lineage_ranks"].fillna("").astype(str)
    #lin2_df["taxid"] = pd.to_numeric(lin2_df["taxid"], errors="coerce").fillna(0).astype(int)
    #lin2_df["full_lineage_taxids"] = lin2_df["full_lineage_taxids"].fillna("")
    #taxid_to_lineage = {int(k): str(v) for k, v in zip(lin2_df["taxid"], lin2_df["full_lineage_taxids"])}
    #taxid_to_lineage_names = dict(zip(
     #   lin2_df["taxid"],
     #   lin2_df["full_lineage_names"]
    #))
    #taxid_to_lineage_ranks = dict(zip(
    #    lin2_df["taxid"],
    #    lin2_df["full_lineage_ranks"]
    #))
    #br["taxon_id"] = pd.to_numeric(br["taxon_id"], errors="coerce").fillna(0).astype(int)
    taxid_to_lineage_names, taxid_to_lineage_ranks, taxid_to_lineage = \
    build_taxonomy_maps(unique_taxids_br, tk_db_dir)
    
    br["full_lineage"] = br["taxon_id"].apply(lambda x: taxid_to_lineage_names.get(int(x), ""))
    br["full_lineage_ranks"] = br["taxon_id"].apply(
        lambda x: taxid_to_lineage_ranks.get(int(x), "")
    )
    br["resolution_level"] = br.apply(classify_taxonomic_resolution, axis=1)
    br["broad_category"] = br.apply(classify_broad_category, axis=1)

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
        "broad_category": "unclassified",
        "reads": unclassified_reads,
        "term_filter": True
    }
    # Append to the bottom of the table
    br = pd.concat([br, pd.DataFrame([unclassified_row])], ignore_index=True)
    br["pc_reads"] = (br["reads"] / filtered_read_counts) * 100
    br["cov_filter"] = br["pc_reads"].astype(float) >= 0.001
    #br_subset = br[['taxon_name', 'taxon_id', 'full_lineage', 'full_lineage_ranks', 'broad_category', 'reads', 'pc_reads', 'term_filter', 'cov_filter']]
    
    #Ensure reads is numeric, then filter detections with >= 10 reads
    br["reads"] = pd.to_numeric(br["reads"], errors="coerce").fillna(0).astype(int)

    br_subset = br[['taxon_name', 'taxon_id', 'full_lineage', 'full_lineage_ranks',
                    'broad_category', 'resolution_level', 'reads', 'pc_reads', 'term_filter', 'cov_filter']]

    # Apply min-reads filter ONLY to non-viral taxa
    # Keep all viral/unclassified-virus rows regardless of reads,
    # but require reads >= 10 for everything else.
    br_subset = filter_viral_rows(br_subset, 0.001)
    #cat = br_subset["broad_category"].astype(str).str.strip().str.lower()

    #viral_labels = {
    #    "viral",
    #    "unclassified virus",
    #    "unclassified viral"
    #}

    #is_viral_taxon = cat.isin(viral_labels)
    #br_subset = br_subset[is_viral_taxon | (br_subset["pc_reads"].astype(float) >= 0.001)].copy()
    out_kraken = os.path.join(args.outdir, f"{sample_name}_kraken_summary.txt")
    br_subset.to_csv(out_kraken, sep="\t", index=False)

if __name__ == "__main__":
    main()

#singularity exec  -B /work/daff_viral_rnaseq/rnaspades/BR_Guar_18 docker://quay.io/biocontainers/pytaxonkit:0.9.1--pyhdfd78af_1 /work/daff_viral_rnaseq/rnaspades/BR_Guar_18/filter_blast.py --blastn_results BR_Guar_18_blastn_header.txt --sample_name BR_Guar_18 --mode ncbi --taxonkit_database_dir ~/.taxonkit
#--sample_name Dataset_1 --bracken /work/Dataset_1_bracken_report.txt --kaiju /work/Dataset_1_kaiju_summary.tsv --taxonkit_database_dir ~/.taxonkit --stats /work/Dataset_1_bbsplit_stats.txt --filter /work/bin/filterKeyWords.txt