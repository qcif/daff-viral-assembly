#!/usr/bin/env python
import argparse
import pandas as pd
import os.path
from Bio import SeqIO


def main():
    ################################################################################
    parser = argparse.ArgumentParser(description="Load blast and coverage stats summary")
    parser.add_argument("--fasta", type=str, required=True, help='provide fasta file')
    parser.add_argument("--sample_name", type=str, required=True, help='provide sample name')
    parser.add_argument("--hmmscan", type=str, required=True, help='provide hmmscan results')
    parser.add_argument("--genomad", type=str, required=True, help='provide genomad results')
    parser.add_argument("--blast", type=str, required=False, help='provide output file name')
    args = parser.parse_args()
    fasta_file = args.fasta
    sample_name = args.sample_name
    genomad = args.genomad
    hmmscan = args.hmmscan
    blast = args.blast
    
    if os.path.getsize(fasta_file) > 0:
        fasta_df = fasta_to_dataframe(fasta_file)
    else:
        fasta_df = pd.DataFrame(columns=["seq_name", "contig_seq"])

    if os.path.getsize(genomad) > 0:
        
        genomad_results = pd.read_csv(genomad, sep="\t", header=0)
        merged_df = pd.merge(fasta_df, genomad_results, on = ['seq_name'], how = 'outer')


    if os.path.getsize(hmmscan) > 0:
        
        hmmscan_results = pd.read_csv(hmmscan, sep="\t", header=0)
        merged2_df = pd.merge(merged_df, hmmscan_results, 
                              left_on = ['seq_name'], 
                              right_on = ['query_name'],
                              how = 'outer')
        
    if os.path.getsize(blast) > 0:
        
        blast_results = pd.read_csv(blast, sep="\t", header=0)
        merged3_df = pd.merge(merged2_df, blast_results, 
                              left_on = ['seq_name'], 
                              right_on = ['qseqid'],
                              how = 'outer')

    num_cols = ["length", "virus_score", "PFAM_total"]
    merged3_df[num_cols] = merged3_df[num_cols].apply(pd.to_numeric, errors="coerce")
    rows_to_check = ["PFAM_total", "virus_score"]

    merged3_df_filt = merged3_df[
        ~(merged3_df[rows_to_check].fillna(0).eq(0).all(axis=1))
    ]
    
    #merged3_df_filt = merged3_df_filt[merged3_df_filt["taxonomy"].str.lower() != "unclassified"]
    merged3_df_filt = merged3_df_filt[merged3_df_filt["length"] >= 500]
    merged3_df_filt = merged3_df_filt[merged3_df_filt["virus_score"] >= 0.9]
    #Keep only contigs with no blast hits.
    merged3_df_filt = merged3_df_filt[
        merged3_df_filt["sacc"].isna() | (merged3_df_filt["sacc"].str.strip() == "")
    ]


    merged3_df_filt["ORFs"] = merged3_df_filt["ORFs"].fillna("NA")

    # Convert RdRp to numeric (True=1, False=0)
    merged3_df_filt["RdRp"] = (
        merged3_df_filt["RdRp"]
        .fillna(False)
        .astype(bool)
        .astype(int)
    )

    # FORCE correct dtypes before sorting
    merged3_df_filt["virus_score"] = pd.to_numeric(
        merged3_df_filt["virus_score"], errors="coerce"
    ).fillna(0)

    merged3_df_filt["PFAM_total"] = pd.to_numeric(
        merged3_df_filt["PFAM_total"], errors="coerce"
    ).fillna(0)


    final_columns_filt = ["seq_name", "contig_seq", "length", "n_genes", "genetic_code", "virus_score", "n_hallmarks", "marker_enrichment", "taxonomy", "ORFs", "RdRp", "PFAM_total"]

    merged3_df_filt = merged3_df_filt[final_columns_filt]
    merged3_df_filt.sort_values(
        by=["virus_score", "RdRp", "PFAM_total"],
        ascending=[False, False, False],
        inplace=True
    )

    output_file = f"{sample_name}_novel_virus_candidates.tsv"
    merged3_df_filt.to_csv(output_file, index=False, sep="\t")

  
# Function to convert FASTA file to DataFrame
def fasta_to_dataframe(fasta_file):
    records = SeqIO.parse(fasta_file, "fasta")
    
    # List to hold the sequence data
    data = []
    
    for record in records:
        # Append ID and sequence to the list
        data.append([record.id, str(record.seq)])
    
    # Create a DataFrame
    df = pd.DataFrame(data, columns=["seq_name", "contig_seq"])
    return df


if __name__ == "__main__":
    main()