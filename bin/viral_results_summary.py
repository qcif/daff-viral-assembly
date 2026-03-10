#!/usr/bin/env python
import argparse
import pandas as pd
import re



# Define PFAMs of interest
PFAMS = [
    "PF00863.26","PF13608.11","PF00998.29","PF00680.26","PF00078.33","PF00767.23",
    "PF08440.16","PF07184.16","PF05533.19","PF01660.23","PF01443.25","PF17646.7",
    "PF01785.23","PF05520.16","PF06922.16","PF03225.20","PF05733.17","PF11479.14",
    "PF11757.14","PF04808.17","PF19223.6","PF00978.27","PF05379.18","PF01307.23",
    "PF20896.3","PF00721.28","PF01107.24","PF01577.23","PF07652.21","PF12524.14", 
    "PF02123.22", "PF00894.25", "PF00729.25", "PF02122.22"
]

# PFAM renaming dictionary
PFAM_RENAME = {
    "PF00863.26": "PF00863.26 (Peptidase_S8_S53)",
    "PF13608.11": "PF13608.11 (Protein P3 of Potyviral polyprotein)",
    "PF00998.29": "PF00998.29 (Viral RNA dependent RNA polymerase)",
    "PF00680.26": "PF00680.26 (Viral RNA-dependent RNA polymerase)",
    "PF00078.33": "PF00078.33 (Reverse transcriptase (RNA-dependent DNA polymerase))",
    "PF00767.23": "PF00767.23 (Potyvirus coat protein)",
    "PF08440.16": "PF08440.16 (Potyviridae polyprotein)",
    "PF07184.16": "PF07184.16 (Citrus tristeza virus P33 protein)",
    "PF05533.19": "PF05533.19 (Peptidase_C42)",
    "PF01660.23": "PF01660.23 (Viral methyltransferase)",
    "PF01443.25": "PF01443.25 (Viral superfamily 1 RNA helicase core domain)",
    "PF17646.7": "PF17646.7 (Closterovirus 1a polyprotein central region)",
    "PF01785.23": "PF01785.23 (Closterovirus coat protein)",
    "PF05520.16": "PF05520.16 (Citrus_P18)",
    "PF06922.16": "PF06922.16 (Citrus tristeza virus P13 protein)",
    "PF03225.20": "PF03225.20 (Viral heat shock protein Hsp90 homologue)",
    "PF05733.17": "PF05733.17 (Tenuivirus/Phlebovirus nucleocapsid protein)",
    "PF11479.14": "PF11479.14 (RNA silencing suppressor P21 C-terminal domain)",
    "PF11757.14": "PF11757.14 (Suppressor of RNA silencing P21-like N-terminal domain)",
    "PF04808.17": "PF04808.17 (Citrus tristeza virus (CTV) P23 protein)",
    "PF19223.6": "PF19223.6 (Chroparavirus methyltransferase)",
    "PF00978.27": "PF00978.27 (RNA dependent RNA polymerase)",
    "PF05379.18": "PF05379.18 (Carlavirus endopeptidase)",
    "PF01307.23": "PF01307.23 (Plant viral movement protein)",
    "PF20896.3": "PF20896.3 (Tomato mosaic virus helicase, N-terminal domain)",
    "PF00721.28": "PF00721.28 (Virus coat protein (TMV like))",
    "PF01107.24": "PF01107.24 (Viral movement protein (MP))",
    "PF01577.23": "PF01577.23 (Potyvirus P1 protease)",
    "PF07652.21": "PF07652.21 (Flavivirus DEAD domain)",
    "PF12524.14": "PF12524.14 (sDNA virus glycoprotein L C terminal)",
    "PF02123.22": "PF02123.22 (Viral RNA-directed RNA-polymerase)",
    "PF00894.25": "PF00894.25 (Luteovirus coat protein)",
    "PF00729.25": "PF00729.25 (Viral coat protein (S domain))",
    "PF02122.22": "PF02122.22 (Peptidase S39)"
}

RDRP_PFAMS = ["PF00998.29","PF00680.26","PF00078.33","PF00978.27","PF02123.22"]

def parse_hmmscan_per_target(path):
    rows = []
    with open(path, "r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = re.split(r"\s+", line.strip(), maxsplit=18)
            rows.append(parts)

    columns = [
        "target_name","target_accession","query_name","query_accession",
        "full_sequence_evalue","full_sequence_score","full_sequence_bias",
        "best_domain_evalue","best_domain_score","best_domain_bias",
        "domain_n_exp","domain_n_reg","domain_n_clu","domain_n_ov",
        "domain_n_env","domain_n_dom","domain_n_rep","domain_n_inc",
        "description_of_target"
    ]

    df = pd.DataFrame(rows, columns=columns)
    df["full_sequence_score"] = df["full_sequence_score"].astype(float)
    return df

def filter_hmmscan(df):
    df_filtered = df[df["full_sequence_score"] >= 50.0]
    df_filtered = df_filtered[df_filtered["target_accession"].isin(PFAMS)]
    df_filtered = df_filtered.drop_duplicates(subset=["query_name","target_accession"])
    return df_filtered

def summarize_domains(df_filtered):
    # Extract ORF number and NODE name
    df_filtered["ORFs"] = df_filtered["query_name"].str.extract(r"_ORF\.(\d+)$")
    df_filtered["query_name"] = df_filtered["query_name"].str.replace(r"_ORF\.\d+$", "", regex=True)

    # Determine which PFAM accessions are present
    pfams_present = sorted(df_filtered["target_accession"].unique())

    # Create count columns dynamically
    for pf in pfams_present:
        df_filtered[pf] = df_filtered["target_accession"].apply(lambda x: 1 if x == pf else 0)

    # Build aggregation dictionary
    agg_dict = {"ORFs": lambda x: ",".join(sorted(x.dropna(), key=lambda v: int(v)))}
    for pf in pfams_present:
        agg_dict[pf] = "sum"

    # Group by NODE
    summary = df_filtered.groupby("query_name").agg(agg_dict).reset_index()

    # Add Total column (sum of PFAM counts)
    summary["PFAM_total"] = summary[pfams_present].sum(axis=1)

    # Remove rows where PFAM count is zero
    summary = summary[summary["PFAM_total"] > 0]

    # If no domains detected → return empty table (no totals row)
    if summary.empty:
        base_cols = ["query_name", "ORFs", "PFAM_total", "RdRp"]
        return pd.DataFrame(columns=base_cols)
    
    # Add RdRp column
    existing_rdRp_pfams = [pf for pf in RDRP_PFAMS if pf in summary.columns]
    #summary["RdRp"] = summary[existing_rdRp_pfams].sum(axis=1) > 0
    summary["RdRp"] = summary[existing_rdRp_pfams].sum(axis=1) > 0 if existing_rdRp_pfams else False


    # Sort by highest number of PFAM occurrences
    summary = summary.sort_values(by="PFAM_total", ascending=False)

    # Rename PFAMs if mapping exists
    rename_dict = {pf: PFAM_RENAME[pf] for pf in pfams_present if pf in PFAM_RENAME}
    summary.rename(columns=rename_dict, inplace=True)

    # Build total summary column
    pfam_cols_final = [rename_dict.get(pf, pf) for pf in pfams_present]
    total_row = summary[pfam_cols_final].sum()
    total_row["PFAM_total"] = summary["PFAM_total"].sum()
    #fill the column row with adequate information
    #total_row["query_name"] = "Total_counts"
    #total_row["ORFs"] = ""
    total_row["Total_counts"] = ""
    #summary = pd.concat([summary, total_row.to_frame().T], ignore_index=True)

    #add column Total
    summary = pd.concat(
        [summary, total_row.to_frame().T],
        ignore_index=True
    )
    return summary
def add_group_max(df, group_col, target_col, new_col):
    df[new_col] = df.groupby(group_col)[target_col].transform("max")
    return df

VIROID_RENAME = {
    "Citrus dwarfing viroid": "Apscaviroid nanocitri",
    "Citrus viroid IIIa": "Apscaviroid nanocitri"
}

def standardize_viroid_names(df, column, df_name=""):
    """Ensure consistent naming for two viroids across any dataframe.

    Converts both known variants to the canonical name and warns if any
    old names remain afterwards.
    """
    if column not in df.columns:
        return df

    df[column] = df[column].replace(VIROID_RENAME)

    remaining = df[column].isin(VIROID_RENAME.keys())
    if remaining.any():
        remaining_vals = df.loc[remaining, column].unique().tolist()
        print(f"WARNING: {df_name} still contains outdated viroid names in '{column}': {remaining_vals}")

    return df


def main():
    ################################################################################
    parser = argparse.ArgumentParser(description="Load blast and coverage stats summary")
    parser.add_argument("--sample_name", type=str, required=True, help='provide sample name')
    parser.add_argument("--blast", type=str, required=True, help='provide blast top hits')
    parser.add_argument("--kraken", type=str, required=True, help='provide kraken file')
    parser.add_argument("--kaiju", type=str, required=True, help='provide kaiju file')
    parser.add_argument("--hmmscan", type=str, required=True, help='provide hmmscan file')
    parser.add_argument("--map2ref", type=str, required=True, help='provide coverage stats to reference')
    args = parser.parse_args()
    sample_name = args.sample_name
    blast = args.blast
    kraken = args.kraken
    kaiju = args.kaiju
    hmmscan = args.hmmscan
    map2ref = args.map2ref
    
    df = parse_hmmscan_per_target(hmmscan)
    df_filtered = filter_hmmscan(df)
    print("Number of matching hits:", len(df_filtered))

    hmm_df = summarize_domains(df_filtered)
    hmm_df.to_csv(sample_name + "_hmm_domain_summary_counts.tsv", sep="\t", index=False)


    blast_df = pd.read_csv(blast, sep="\t", dtype=str)
    kraken_df = pd.read_csv(kraken, sep="\t", dtype=str)
    kaiju_df = pd.read_csv(kaiju, sep="\t", dtype=str)
    map2ref_df = pd.read_csv(map2ref, sep="\t", dtype=str)
    print(blast_df.head())

    #filter blast results to only those that passed all filters inc best contig per species filter
    blast_df["pident"] = pd.to_numeric(blast_df["pident"], errors="coerce")
    blast_df["qlen"] = pd.to_numeric(blast_df["qlen"], errors="coerce")
    blast_df["mapping_read_count"] = pd.to_numeric(blast_df["mapping_read_count"], errors="coerce")
    blast_df["pc_mapping_reads"] = pd.to_numeric(blast_df["pc_mapping_reads"], errors="coerce")
    blast_df = add_group_max(blast_df, "species", "pident", "max_pident_spp")
    blast_df = add_group_max(blast_df, "species", "qlen", "max_qlen_spp")
    blast_df = add_group_max(blast_df, "species", "mapping_read_count", "max_mapping_read_count_spp")
    blast_df = add_group_max(blast_df, "species", "pc_mapping_reads", "max_pc_mapping_reads_spp")


    #only reports PFAMs for contigs that passed all filters and were the best contig per species
    #chnage to stre whether there was a PFAM hit for any contig per species that passed filters, not just the best contig per species
    blast_df2 = blast_df.merge(
        hmm_df[["query_name", "PFAM_total", "RdRp"]],
        left_on="qseqid",
        right_on="query_name",
        how="left"   # left join keeps all rows in summary_df
    )
    blast_df2 = add_group_max(blast_df2, "species", "PFAM_total", "max_PFAM_total_spp")
    # ensure RdRp is boolean
    blast_df2["RdRp"] = blast_df2["RdRp"].fillna(False)

    # species-level flag: any contig of this species has RdRp
    blast_df2["species_has_RdRp"] = (
        blast_df2.groupby("species")["RdRp"]
        .transform("any")
    )
    

    filtered_blast_df = blast_df2[
        (blast_df2["term_filter"].astype(str) == "True") &
        (blast_df2["cov_filter"].astype(str) == "True") &
        (blast_df2["best_contig_per_sp_filter"].astype(str) == "True")
    ]

    # Standardize species naming across all tables for the same viroid
    filtered_blast_df = standardize_viroid_names(filtered_blast_df, "species", df_name="filtered_blast_df")
    kraken_df = standardize_viroid_names(kraken_df, "taxon_name", df_name="kraken_df")
    kaiju_df = standardize_viroid_names(kaiju_df, "taxon_name", df_name="kaiju_df")
    map2ref_df = standardize_viroid_names(map2ref_df, "taxon_name", df_name="map2ref_df")

    #print(kraken_df.head())
    #print(kaiju_df.head())
    summary_df = pd.DataFrame()
    
    
    # Mask Kraken/Kaiju species to only viral entries and filter for those that passed the term filter
    kraken_viral = kraken_df.loc[
        (kraken_df["broad_categories"] == "viral") & (kraken_df["term_filter"].str.lower() == "true"),
        "taxon_name"
    ]
    kaiju_viral  = kaiju_df.loc[
        (kaiju_df["broad_categories"] == "viral") & (kraken_df["term_filter"].str.lower() == "true"),
        "taxon_name"
    ]

    # Combine vertically and get unique
    summary_df["taxon"] = pd.Series(
        pd.unique(
            pd.concat([
                filtered_blast_df["species"],
                kraken_viral,
                kaiju_viral
            ], ignore_index=True)
        )
    )
    # Add boolean column: True if species appeared in blast
    blast_species_set = set(filtered_blast_df["species"])
    kraken_species_set = set(kraken_df["taxon_name"])
    kaiju_species_set = set(kaiju_df["taxon_name"])
    summary_df["taxon_in_blast"] = summary_df["taxon"].isin(blast_species_set)
    summary_df["taxon_in_kraken"] = summary_df["taxon"].isin(kraken_species_set)
    summary_df["taxon_in_kaiju"] = summary_df["taxon"].isin(kaiju_species_set)
    #print(summary_df)

    merged_df = summary_df.merge(
        filtered_blast_df[["species", "qseqid", "qlen", "sacc", "pident", "bitscore", "evalue", "contig_seq", "ncontigs_per_spp", 
                           "max_pident_spp","max_qlen_spp", "max_mapping_read_count_spp",  "max_pc_mapping_reads_spp", "total_score_spp", "mapping_read_count", "pc_mapping_reads", "mean_depth", "pc_cov_30X", 
                            "mean_mapping_quality", "read_count_flag", "mean_depth_flag", "30x_cov_flag", "mean_mq_flag", "species_has_RdRp", "max_PFAM_total_spp",
                            "total_conf_score","normalised_conf_score"]],
        left_on="taxon",
        right_on="species",
        how="left"   # left join keeps all rows in summary_df
    )
    merged_df.rename(columns={"qlen": "contig_length",
                               "contig_seq": "contig_seq"}, inplace=True)


    merged_df2 = merged_df.merge(
        kaiju_df[["taxon_name", "reads", "pc_reads"]],
        left_on="taxon",
        right_on="taxon_name",
        how="left"   # left join keeps all rows in summary_df
    )
    merged_df2[["reads", "pc_reads"]] = (
        merged_df2[["reads", "pc_reads"]].fillna(0)
    )
    # Rename the columns
    merged_df2.rename(columns={"pc_reads": "kaiju_pc_reads",
                               "reads": "kaiju_reads"}, inplace=True)


    merged_df3 = merged_df2.merge(
        kraken_df[["taxon_name", "reads", "pc_reads"]],
        left_on="taxon",
        right_on="taxon_name",
        how="left"   # left join keeps all rows in summary_df
    )
    merged_df3[["reads", "pc_reads"]] = (
        merged_df3[["reads", "pc_reads"]].fillna(0)
    )
    # Rename the columns
    merged_df3.rename(columns={"pc_reads": "kraken_pc_reads",
                               "reads": "kraken_reads"}, inplace=True)

   
    map2ref_df["ref_count"] = map2ref_df.groupby(["taxon_name", "sacc"])["sacc"].transform("count")
    map2ref_df = add_group_max(map2ref_df, "taxon_name", "pc_cov_30X", "max_pc_cov_30X_spp")
    map2ref_df = add_group_max(map2ref_df, "taxon_name", "normalised_conf_score", "max_normalised_conf_score_spp")
    map2ref_df = add_group_max(map2ref_df, "taxon_name", "reference_length", "max_reference_length_spp")
    map2ref_df_unique = (
        map2ref_df[
            ["taxon_name",
            "max_reference_length_spp",
            "max_normalised_conf_score_spp",
            "max_pc_cov_30X_spp"]
            ]
            .drop_duplicates(subset=["taxon_name"])
        )
    
    merged_df4 = merged_df3.merge(
        map2ref_df_unique[["taxon_name", "max_reference_length_spp", "max_normalised_conf_score_spp", "max_pc_cov_30X_spp"]],
        left_on="taxon",
        right_on="taxon_name",
        how="left"   # left join keeps all rows in summary_df
    )

    final_columns_filt = ["taxon","kraken_reads", "kraken_pc_reads", 
                          "max_mapping_read_count_spp",  "max_pc_mapping_reads_spp",
                          "kaiju_reads", "kaiju_pc_reads", 
                          "mapping_read_count","pc_mapping_reads",
                          "ncontigs_per_spp", "qseqid", "contig_seq", "max_pident_spp", "max_qlen_spp",
                          "mean_depth", "pc_cov_30X", "normalised_conf_score", "max_PFAM_total_spp", "species_has_RdRp", 
                          "max_reference_length_spp", "max_normalised_conf_score_spp", "max_pc_cov_30X_spp"]

    merged_df4 = merged_df4[final_columns_filt]
    text_cols = ["qseqid", "contig_seq","species_has_RdRp" ]
    #text_cols = ["qseqid", "contig_seq", "sacc", "read_count_flag", "mean_depth_flag", "30x_cov_flag", "mean_mq_flag", "RdRp" ]
    num_cols = ["ncontigs_per_spp","mapping_read_count", "pc_mapping_reads", "mean_depth", "pc_cov_30X", "normalised_conf_score" , "max_PFAM_total_spp",  "max_pident_spp", "max_qlen_spp"]

    # fill values
    merged_df4 = merged_df4[merged_df4["contig_seq"].notna()]
    merged_df4[text_cols] = merged_df4[text_cols].fillna("NA")
    merged_df4[num_cols] = merged_df4[num_cols].fillna(0)
    
    merged_df4[num_cols] = merged_df4[num_cols].apply(pd.to_numeric, errors="coerce")

    merged_df4.sort_values(
        by=["mapping_read_count", "normalised_conf_score"],
        ascending=[False, False],
        inplace=True
    )
    
    output_file = f"{sample_name}_summary_viral_results.tsv"
    merged_df4.to_csv(output_file, index=False, sep="\t")

if __name__ == "__main__":
    main()