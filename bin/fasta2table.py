#!/usr/bin/env python
import argparse
import pandas as pd
import os.path
from Bio import SeqIO


def main():
    ################################################################################
    parser = argparse.ArgumentParser(description="Load blast and coverage stats summary")
    parser.add_argument("--fasta", type=str, required=True, help='provide fasta file')
    parser.add_argument("--sample", type=str, required=True, help='provide sample name')
    parser.add_argument("--tophits", type=str, required=True, help='provide blast top hits')
    args = parser.parse_args()
    fasta_file = args.fasta
    sample_name = args.sample
    blast = args.tophits
    if os.path.getsize(blast) > 0:
        fasta_df = fasta_to_dataframe(fasta_file)
    else:
        fasta_df = pd.DataFrame(columns=["qseqid", "contig_seq"])

    if os.path.getsize(blast) > 0:
        blastn_results = pd.read_csv(blast, sep="\t", header=0)
        blastn_results.drop(['sample_name'], axis=1, inplace=True)
        merged_df = pd.merge(fasta_df, blastn_results, on = ['qseqid'], how = 'inner')
        merged_df.insert(0, "sample_name", sample_name)
        #merged_df['n_read_cont_cluster'] = merged_df['qseqid'].str.split('_').str[5].astype(int)
       #merged_df['n_read_cont_cluster'] = merged_df['n_read_cont_cluster'].str.replace("RC","").astype(int)
        merged_df = merged_df.sort_values(["best_contig_per_sp_filter"], ascending=[False])

        # merged_df.to_csv(str(sample_name) + "_blastn_top_hits.txt", index=None, sep="\t")
        merged_df.to_csv(os.path.basename(blast).replace("_top_viral_hits.txt", "_top_viral_hits_with_contigs.txt"), index=None, sep="\t")
        filtered_df = merged_df[merged_df["term_filter"] & merged_df["cov_filter"]]
        filtered_df.to_csv(os.path.basename(blast).replace("_top_viral_hits.txt", "_top_viral_hits_filtered_with_contigs.txt"), index=None, sep="\t")
        
        # Extract column 4
        col4 = filtered_df.iloc[:, 3].astype(str)  # ensure strings

        # Replace spaces with underscores
        col4 = col4.str.replace(" ", "_")

        # Flatten values and get unique, sorted IDs
        unique_ids = list(dict.fromkeys(col4))

        # Save to file, one ID per line
        with open(sample_name + "_ids_to_retrieve.txt", "w") as f:
            for uid in unique_ids:
                f.write(f"{uid}\n")

    else:
        print("DataFrame is empty!")
        for col in ["sample_name", "qseqid", "sacc", "alignment_length", "evalue", "bitscore", "pident", "mismatch",
                     "gapopen", "qstart", "qend", "qlen", "sstart", "send", "slen", "sstrand",
                     "qcovhsp", "staxids", "qseq", "sseq", "qcovs", 
                     "species_updated", "RNA_type", "stitle", "full_lineage", "ncontigs", "total_score", "term_filter", "cov_filter", "best_contig_per_sp_filter"]:
            if col not in fasta_df.columns:
                fasta_df[col] = None
                fasta_df.to_csv(os.path.basename(blast).replace("_top_viral_hits.txt", "_top_viral_hits_with_contigs.txt"), index=None, sep="\t")
                filtered_df[col] = None
                filtered_df.to_csv(os.path.basename(blast).replace("_top_viral_hits.txt", "_top_viral_hits_filtered_with_contigs.txt"), index=None, sep="\t")
# Function to convert FASTA file to DataFrame
def fasta_to_dataframe(fasta_file):
    records = SeqIO.parse(fasta_file, "fasta")
    
    # List to hold the sequence data
    data = []
    
    for record in records:
        # Append ID and sequence to the list
        data.append([record.id, str(record.seq)])
    
    # Create a DataFrame
    df = pd.DataFrame(data, columns=["qseqid", "contig_seq"])
    return df


if __name__ == "__main__":
    main()