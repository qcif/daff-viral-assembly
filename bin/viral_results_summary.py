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
    
    # Add RdRp column
    existing_rdRp_pfams = [pf for pf in RDRP_PFAMS if pf in summary.columns]
    summary["RdRp"] = summary[existing_rdRp_pfams].sum(axis=1) > 0


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


def main():
    ################################################################################
    parser = argparse.ArgumentParser(description="Load blast and coverage stats summary")

    parser.add_argument("--sample_name", type=str, required=True, help='provide sample name')
    parser.add_argument("--blast", type=str, required=True, help='provide blast top hits')
    parser.add_argument("--kraken", type=str, required=True, help='provide fasta file')
    parser.add_argument("--kaiju", type=str, required=True, help='provide fasta file')
    parser.add_argument("--hmmscan", type=str, required=True, help='provide hmmscan file')
    #parser.add_argument("--map2ref", type=str, required=True, help='provide coverage stats to reference')
    args = parser.parse_args()
    sample_name = args.sample_name
    blast = args.blast
    kraken = args.kraken
    kaiju = args.kaiju
    hmmscan = args.hmmscan
    #map2ref = args.map2ref
    
    df = parse_hmmscan_per_target(hmmscan)
    df_filtered = filter_hmmscan(df)
    print("Number of matching hits:", len(df_filtered))

    hmm_df = summarize_domains(df_filtered)
    hmm_df.to_csv(sample_name + "_hmm_domain_summary_counts.tsv", sep="\t", index=False)


    blast_df = pd.read_csv(blast, sep="\t", dtype=str)
    kraken_df = pd.read_csv(kraken, sep="\t", dtype=str)
    kaiju_df = pd.read_csv(kaiju, sep="\t", dtype=str)
    #map2ref_df = pd.read_csv(map2ref, sep="\t", dtype=str)

    print(blast_df.head())

    
    filtered_blast_df = blast_df[
        (blast_df["term_filter"].astype(str) == "True") &
        (blast_df["cov_filter"].astype(str) == "True") &
        (blast_df["best_contig_per_sp_filter"].astype(str) == "True")
    ]


    # Standardize species naming
    filtered_blast_df["species"] = filtered_blast_df["species"].replace(
        {"Citrus viroid IIIa": "Apscaviroid nanocitri"}
    )

    kraken_df["taxon_name"] = kraken_df["taxon_name"].replace(
        {"Citrus dwarfing viroid": "Apscaviroid nanocitri"}
    )

    kaiju_df["taxon_name"] = kaiju_df["taxon_name"].replace(
        {"Citrus dwarfing viroid": "Apscaviroid nanocitri"}
    )

    #print(kraken_df.head())
    #print(kaiju_df.head())
    summary_df = pd.DataFrame()
    
    # Mask Kraken/Kaiju species to only viral entries
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
        filtered_blast_df[["species", "qseqid", "qlen", "sacc", "pident", "bitscore", "evalue", "contig_seq", "ncontigs_per_spp", "total_score_spp", "mapping_read_count", "pc_mapping_reads", "mean_depth", "pc_cov_30X",  
        "mean_mapping_quality", "read_count_flag", "mean_depth_flag", "30x_cov_flag", "mean_mq_flag", "total_conf_score","normalised_conf_score"]],
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
    
    

    merged_df4 = merged_df3.merge(
        hmm_df[["query_name", "PFAM_total", "RdRp"]],
        left_on="qseqid",
        right_on="query_name",
        how="left"   # left join keeps all rows in summary_df
    )

    final_columns_filt = ["taxon", "taxon_in_blast", "taxon_in_kraken", "taxon_in_kaiju", "kraken_reads", "kraken_pc_reads", "kaiju_reads", "kaiju_pc_reads", "ncontigs_per_spp", "qseqid", "contig_seq", "contig_length", "pident", "bitscore", "evalue", "sacc", "mapping_read_count","pc_mapping_reads", "mean_depth", "pc_cov_30X",  
        "mean_mapping_quality", "read_count_flag", "mean_depth_flag", "30x_cov_flag", "mean_mq_flag",  "total_conf_score","normalised_conf_score", "PFAM_total", "RdRp"]

    merged_df4 = merged_df4[final_columns_filt]
    
    text_cols = ["qseqid", "contig_seq", "sacc", "read_count_flag", "mean_depth_flag", "30x_cov_flag", "mean_mq_flag", "RdRp" ]
    num_cols = ["ncontigs_per_spp", "contig_length", "pident", "bitscore", "evalue", "mapping_read_count", "pc_mapping_reads", "mean_depth", "pc_cov_30X", "mean_mapping_quality", "total_conf_score", "normalised_conf_score" , "PFAM_total"]

    # fill values
    merged_df4_filt = merged_df4[merged_df4["contig_seq"].notna()]
    merged_df4_filt[text_cols] = merged_df4_filt[text_cols].fillna("NA")
    merged_df4_filt[num_cols] = merged_df4_filt[num_cols].fillna(0)
    
    merged_df4_filt[num_cols] = merged_df4_filt[num_cols].apply(pd.to_numeric, errors="coerce")

    merged_df4_filt.sort_values(
        by=["mapping_read_count", "normalised_conf_score"],
        ascending=[False, False],
        inplace=True
    )
    
    output_file = f"{sample_name}_summary_viral_results.tsv"
    merged_df4_filt.to_csv(output_file, index=False, sep="\t")
    
    #merged_df4 = merged_df3.merge(
    #    map2ref_df[["taxon_name", "consensus_seq", "sacc", "reference_length", "reads", "pc_reads", "mean_depth", "pc_cov_30X", "mean_mapping_quality", "normalised_conf_score" ]],
    #    left_on=["taxon","sacc"],
    #    right_on=["taxon_name","sacc"],
    #    how="left"   # left join keeps all rows in summary_df
    #)

    #merged_df4.sort_values(
    #    by=["contig_length", "kraken_pc_reads", "kaiju_pc_reads"], 
    #    ascending=[False, False, False], 
    #    inplace=True
    #)
    # Rename the columns
    #merged_df4.rename(columns={"consensus_seq": "ref_consensus_seq"}, inplace=True)
    #merged_df4.fillna({
    #    "reference_length": 0,
    #    "reads": 0,
    #    "species": "NA",
    #    "genus": "NA",
    #    "description": "NA"
    #}, inplace=True)

    #final_columns_filt = ["taxon", "taxon_in_blast", "taxon_in_kraken", "taxon_in_kaiju", "kraken_reads", "kraken_pc_reads", "kaiju_reads", "kaiju_pc_reads", "ncontigs_per_spp", "qseqid", "contig_seq", "contig_length", "pident", "bitscore", "evalue", "sacc", "reference_length", "ref_consensus_seq", "reads", "pc_reads", "mean_depth", "pc_cov_30X", "mean_mapping_quality", "normalised_conf_score"]
    
    #merged_df4 = merged_df4[final_columns_filt]
    #output_file = f"{sample_name}_summary_viral_results.tsv"
    #merged_df4.to_csv(output_file, index=False, sep="\t")


if __name__ == "__main__":
    main()