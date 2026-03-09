#!/usr/bin/env python
"""
Convert FASTA sequences and BLAST results to annotated tables.
Supports two modes: viral hits and reference mapping.
"""
import argparse
import pandas as pd
import os.path
import sys
from Bio import SeqIO

def fasta_to_dataframe(fasta_file):
    """Convert FASTA file to DataFrame with columns [qseqid, contig_seq]."""
    records = SeqIO.parse(fasta_file, "fasta")
    # List to hold the sequence data
    data = []
    
    try:
        for record in records:
             # Append ID and sequence to the list
            data.append([record.id, str(record.seq)])
    except Exception as e:
        print(f"Error parsing FASTA file {fasta_file}: {e}", file=sys.stderr)
        sys.exit(1)
    
    if not data:
        print(f"Warning: No records found in FASTA file {fasta_file}", file=sys.stderr)
    
    return pd.DataFrame(data, columns=["qseqid", "contig_seq"])


def process_denovo_contigs(fasta_df, blastn_results, sample_name, blast_file):
    """Process BLAST results for viral hits."""
    try:
        # Validate expected columns
        required_cols = ["qseqid", "term_filter", "cov_filter", "best_contig_per_sp_filter"]
        missing_cols = [col for col in required_cols if col not in blastn_results.columns]
        if missing_cols:
            raise ValueError(f"Missing required columns for viral mode: {missing_cols}")
        
        # Remove sample_name column if it exists
        if "sample_name" in blastn_results.columns:
            blastn_results = blastn_results.drop(['sample_name'], axis=1)
        
        # Merge FASTA with BLAST results
        merged_df = pd.merge(fasta_df, blastn_results, on=['qseqid'], how='inner')
        
        if merged_df.empty:
            print(f"Warning: No matching sequences found after merge", file=sys.stderr)
        
        # Add sample name column and sort
        merged_df.insert(0, "sample_name", sample_name)
        merged_df = merged_df.sort_values(["best_contig_per_sp_filter"], ascending=[False])
        
        # Save full results
        output_file = os.path.basename(blast_file).replace("_top_viral_hits.txt", "_top_viral_hits_with_contigs.txt")
        merged_df.to_csv(output_file, index=False, sep="\t")
        
        # Apply filters
        filtered_df = merged_df[merged_df["term_filter"] & merged_df["cov_filter"]]
        
        # Save filtered results
        filtered_output = os.path.basename(blast_file).replace("_top_viral_hits.txt", "_top_viral_hits_filtered_with_contigs.txt")
        filtered_df.to_csv(filtered_output, index=False, sep="\t")
        
        # Extract and save reference IDs
        if not filtered_df.empty and len(filtered_df.columns) > 3:
            col4 = filtered_df.iloc[:, 3].astype(str)
            col4 = col4.str.replace(" ", "_")
            unique_ids = list(dict.fromkeys(col4))
            
            with open(f"{sample_name}_ref_ids_to_retrieve.txt", "w") as f:
                for uid in unique_ids:
                    f.write(f"{uid}\n")
        else:
            # Create empty file if no filtered results
            with open(f"{sample_name}_ref_ids_to_retrieve.txt", "w") as f:
                pass
                
    except Exception as e:
        print(f"Error processing viral hits: {e}", file=sys.stderr)
        sys.exit(1)


def process_reference(fasta_df, blastn_results, sample_name, blast_file):
    """Process BLAST results for reference mapping."""
    try:
        # Validate expected columns
        required_cols = ["sacc", "normalised_conf_score"]
        missing_cols = [col for col in required_cols if col not in blastn_results.columns]
        if missing_cols:
            raise ValueError(f"Missing required columns for reference mode: {missing_cols}")
        
        # Extract accession from qseqid
        fasta_df = fasta_df.copy()
        fasta_df["sacc"] = fasta_df["qseqid"].str.split(".").str[0]
        
        # Merge on accession
        merged_df = pd.merge(fasta_df, blastn_results, on=['sacc'], how='inner')
        
        if merged_df.empty:
            print(f"Warning: No matching sequences found after merge", file=sys.stderr)
        
        # Rename and clean up
        merged_df.insert(0, "sample_name", sample_name)
        merged_df = merged_df.rename(columns={"contig_seq": "consensus_seq"})
        merged_df = merged_df.drop("qseqid", axis=1)
        merged_df = merged_df.sort_values(["normalised_conf_score"], ascending=[False])
        
        # Save results
        output_file = os.path.basename(blast_file).replace("_reference_with_cov_stats.txt", "_reference_with_cov_stats_final.txt")
        merged_df.to_csv(output_file, index=False, sep="\t")
        
    except Exception as e:
        print(f"Error processing reference stats: {e}", file=sys.stderr)
        sys.exit(1)


def create_empty_outputs(sample_name, blast_file, mode="viral"):
    """Create empty output files with proper structure."""
    expected_columns = [
        "sample_name", "qseqid", "sacc", "alignment_length", "evalue", "bitscore",
        "pident", "mismatch", "gapopen", "qstart", "qend", "qlen", "sstart", "send",
        "slen", "sstrand", "qcovhsp", "staxids", "qseq", "sseq", "qcovs", "species",
        "RNA_type", "stitle", "full_lineage", "ncontigs_per_sacc", "ncontigs_per_spp",
        "total_score", "term_filter", "cov_filter", "best_contig_per_sp_filter"
    ]
    
    empty_df = pd.DataFrame(columns=expected_columns)
    
    if mode == "viral":
        output_file = os.path.basename(blast_file).replace("_top_viral_hits.txt", "_top_viral_hits_with_contigs.txt")
        filtered_file = os.path.basename(blast_file).replace("_top_viral_hits.txt", "_top_viral_hits_filtered_with_contigs.txt")
        empty_df.to_csv(output_file, index=False, sep="\t")
        empty_df.to_csv(filtered_file, index=False, sep="\t")
        with open(f"{sample_name}_ref_ids_to_retrieve.txt", "w") as f:
            pass
    elif mode == "reference":
        output_file = os.path.basename(blast_file).replace("_reference_with_cov_stats.txt", "_reference_with_cov_stats_final.txt")
        empty_df.to_csv(output_file, index=False, sep="\t")


def main():
    parser = argparse.ArgumentParser(
        description="Convert FASTA and BLAST results to annotated tables",
        epilog="Examples:\n"
               "  fasta2table.py --fasta contigs.fa --sample sample1 --tophits hits.txt --mode contigs\n"
               "  fasta2table.py --fasta refs.fa --sample sample1 --tophits stats.txt --mode reference"
    )
    
    parser.add_argument("--fasta", type=str, required=True, help='FASTA file path')
    parser.add_argument("--sample", type=str, required=True, help='Sample name')
    parser.add_argument("--tophits", type=str, required=True, help='BLAST results file')
    parser.add_argument("--mode", type=str, required=True, choices=["contigs", "reference"],
                        help='Processing mode: contigs or reference (with coverage stats)')
    args = parser.parse_args()
    
    fasta_file = args.fasta
    sample_name = args.sample
    blast_file = args.tophits
    mode = args.mode
    
    # Load FASTA data
    if os.path.getsize(fasta_file) > 0:
        fasta_df = fasta_to_dataframe(fasta_file)
    else:
        print(f"Warning: FASTA file is empty: {fasta_file}", file=sys.stderr)
        fasta_df = pd.DataFrame(columns=["qseqid", "contig_seq"])
    
    # Load and process BLAST data
    if os.path.getsize(blast_file) > 0:
        try:
            blastn_results = pd.read_csv(blast_file, sep="\t", header=0)
        except Exception as e:
            print(f"Error reading BLAST file {blast_file}: {e}", file=sys.stderr)
            sys.exit(1)
        
        if blastn_results.empty:
            print(f"Warning: BLAST file contains no data rows", file=sys.stderr)
            create_empty_outputs(sample_name, blast_file, mode=mode)
            return
        
        # Route to appropriate handler based on mode
        if mode == "contigs":
            process_denovo_contigs(fasta_df, blastn_results, sample_name, blast_file)
        elif mode == "reference":
            process_reference(fasta_df, blastn_results, sample_name, blast_file)
    else:
        print(f"Warning: BLAST file is empty: {blast_file}", file=sys.stderr)
        create_empty_outputs(sample_name, blast_file, mode=mode)

if __name__ == "__main__":
    main()