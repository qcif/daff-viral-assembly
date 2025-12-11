#!/usr/bin/env python
import argparse
import pandas as pd

def combine_unique(a, b, c):
    return list(pd.unique([a, b, c]))

def main():
    ################################################################################
    parser = argparse.ArgumentParser(description="Load blast and coverage stats summary")

    parser.add_argument("--sample_name", type=str, required=True, help='provide sample name')
    parser.add_argument("--blast", type=str, required=True, help='provide blast top hits')
    parser.add_argument("--kraken", type=str, required=True, help='provide fasta file')
    parser.add_argument("--kaiju", type=str, required=True, help='provide fasta file')
    parser.add_argument("--map2ref", type=str, required=True, help='provide coverage stats to reference')
    args = parser.parse_args()
    sample_name = args.sample_name
    blast = args.blast
    kraken = args.kraken
    kaiju = args.kaiju
    map2ref = args.map2ref
    
    blast_df = pd.read_csv(blast, sep="\t", dtype=str)
    kraken_df = pd.read_csv(kraken, sep="\t", dtype=str)
    kaiju_df = pd.read_csv(kaiju, sep="\t", dtype=str)
    #map2ref_df = pd.read_csv(map2ref, sep="\t", dtype=str)

    print(blast_df.head())
    filtered_blast_df = blast_df[
        (blast_df["term_filter"].astype(str) == "True") &
        #(blast_df["cov_filter"].astype(str) == "True") &
        (blast_df["best_contig_per_sp_filter"].astype(str) == "True")
    ]
    print(kraken_df.head())
    print(kaiju_df.head())
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
    print(summary_df)

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
    final_columns_filt = ["taxon", "taxon_in_blast", "taxon_in_kraken", "taxon_in_kaiju", "kraken_reads", "kraken_pc_reads", "kaiju_reads", "kaiju_pc_reads", "ncontigs_per_spp", "qseqid", "contig_seq", "contig_length", "pident", "bitscore", "evalue", "sacc", "mapping_read_count","pc_mapping_reads", "mean_depth", "pc_cov_30X",  
        "mean_mapping_quality", "read_count_flag", "mean_depth_flag", "30x_cov_flag", "mean_mq_flag", "total_conf_score","normalised_conf_score"]
    merged_df3 = merged_df3[final_columns_filt]
    
    text_cols = ["qseqid", "contig_seq", "sacc", "read_count_flag", "mean_depth_flag", "30x_cov_flag", "mean_mq_flag" ]
    num_cols = ["ncontigs_per_spp", "contig_length", "pident", "bitscore", "evalue", "mapping_read_count", "pc_mapping_reads", "mean_depth", "pc_cov_30X", "mean_mapping_quality", "total_conf_score", "normalised_conf_score" ]

    # fill values
    
    merged_df3[text_cols] = merged_df3[text_cols].fillna("NA")
    merged_df3[num_cols] = merged_df3[num_cols].fillna(0)
    
    merged_df3[num_cols] = merged_df3[num_cols].apply(pd.to_numeric, errors="coerce")
    print(merged_df3.dtypes)
    merged_df3.sort_values(
        by=["mapping_read_count", "normalised_conf_score"],
        ascending=[False, False],
        inplace=True
    )

    output_file = f"{sample_name}_summary_viral_results.tsv"
    merged_df3.to_csv(output_file, index=False, sep="\t")
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