#!/usr/bin/env python
import argparse
import pandas as pd
import numpy as np
from functools import reduce
from Bio import SeqIO
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import collections
import json
import re


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

def parse_args():
    parser = argparse.ArgumentParser(description="Load blast and coverage stats summary")
    parser.add_argument("--sample", type=str, required=True, help='Provide sample name')
    parser.add_argument("--blastn_results", type=str, required=True)
    parser.add_argument("--bbsplit_stats", type=str, required=True)
    parser.add_argument("--bed", type=str, required=True)
    parser.add_argument("--coverage", type=str, required=True)
    parser.add_argument("--mapping_quality", type=str, required=True)
    return parser.parse_args()


def read_filtered_read_count(fastp_path):
    filtered_reads = 0
    with open(fastp_path) as f:
        data = json.load(f)
    filtered_reads = data["summary"]["after_filtering"]["total_reads"]

    return filtered_reads
        

def load_and_prepare_data(coverage_path, bed_path, mapping_quality, filtered_read_counts):
    

    samtools_cov = pd.read_csv(coverage_path, sep="\t", usecols=["#rname", "endpos", "numreads", "meandepth"], header=0)
    samtools_cov.rename(columns={
        "#rname": "qseqid",
        "endpos": "query_match_length",
        "numreads": "qseq_mapping_read_count",
        "meandepth": "qseq_mean_depth"
    }, inplace=True)
    samtools_cov['qseq_pc_mapping_read'] = samtools_cov['qseq_mapping_read_count'] / filtered_read_counts * 100

    mosdepth = pd.read_csv(bed_path, sep="\t", header=0)
    mosdepth.columns = ["qseqid", "start", "end", "region", "base_counts_at_depth_30X"]
    mosdepth['qseq_pc_cov_30X'] = np.where(
        mosdepth['base_counts_at_depth_30X'] > mosdepth['end'],
        100,
        mosdepth['base_counts_at_depth_30X'] / mosdepth['end'] * 100
    )
    mosdepth['qseq_pc_cov_30X'] = mosdepth['qseq_pc_cov_30X'].round(1)

    mq = pd.read_csv(mapping_quality, sep="\t", header=None)
    mq.columns = ["qseqid", "mean_MQ"]
    return samtools_cov, mosdepth[["qseqid", "qseq_pc_cov_30X"]],mq


def merge_dataframes(samtools_cov, mosdepth_df, mq_df):
    return reduce(
        lambda left, right: pd.merge(left, right, on="qseqid", how='outer').fillna(0),
        [samtools_cov, mosdepth_df, mq_df]
    )


def apply_qc_flags(df):
    #Conditions:
    #GREEN: If sgi != 0 and the qseq_pc_cov_30X is >=90
    #ORANGE: If sgi != 0 and the qseq_pc_cov_30X is between 75 and 90.
    #RED: If sgi != 0 and the qseq_pc_cov_30X is <75.
    #GREY: If sgi == 0.

    df['30X_COVERAGE_FLAG'] = np.select(
        [
            (df['qseqid'] != 0) & (df['qseq_pc_cov_30X'] >= 90),
            (df['qseqid'] != 0) & (df['qseq_pc_cov_30X'] >= 75) & (df['qseq_pc_cov_30X'] < 90),
            (df['qseqid'] != 0) & (df['qseq_pc_cov_30X'] < 75),
            (df['qseqid'].isin([0, None, '', '0', '-']))
        ],
        ['GREEN', 'ORANGE', 'RED', 'GREY'],
        default=""
    )

    #######MAPPED_READ_COUNT_FLAG#######
    #Conditions:
    #GREEN: If sgi != 0 and the qseq_mapping_read_count is >=1000
    #ORANGE: If sgi != 0 and the qseq_mapping_read_count is between 200 and 1000.
    #RED: If sgi != 0 and the qseq_mapping_read_count is <200.
    #GREY: If sgi == 0.
    df['MAPPED_READ_COUNT_FLAG'] = np.select(
        [
            (df['qseqid'] != 0) & (df['qseq_mapping_read_count'] >= 1000),
            (df['qseqid'] != 0) & (df['qseq_mapping_read_count'] >= 200) & (df['qseq_mapping_read_count'] < 1000),
            (df['qseqid'] != 0) & (df['qseq_mapping_read_count'] < 200),
            (df['qseqid'].isin([0, None, '', '0', '-']))
        ],
        ['GREEN', 'ORANGE', 'RED', 'GREY'],
        default=""
    )

    # Mean coverage flag
    df['MEAN_COVERAGE_FLAG'] = np.select(
        [
            (df['qseqid'] != 0) & (df['qseq_mean_depth'] >= 500),
            (df['qseqid'] != 0) & (df['qseq_mean_depth'] >= 100) & (df['qseq_mean_depth'] < 500),
            (df['qseqid'] != 0) & (df['qseq_mean_depth'] < 100),
            (df['qseqid'].isin([0, None, '', '0', '-']))
        ],
        ['GREEN', 'ORANGE', 'RED', 'GREY'],
        default=""
    )

    # Mean mapping quality flag
    df['MEAN_MQ_FLAG'] = np.where(
        (df['qseqid'] != 0) & 
        (df['mean_MQ'] >= 30),
        "GREEN",
        np.where((df['qseqid'] != 0) & 
                (df['mean_MQ'] < 30) & 
                (df['mean_MQ'] >= 10),
                "ORANGE",
            np.where((df['qseqid'] != 0) &
                (df['mean_MQ'] < 10),
                "RED",
                np.where(df['qseqid'].isin([0, None, '', '0', '-']),
                    "GREY",
                    ""
                )
            )
        )
    )

    flag_score_map = {
        'GREEN': 2,
        'ORANGE': 1,
        'RED': 0,
        'GREY': 0
    }
    
    flag_columns = [
        '30X_COVERAGE_FLAG',
        'MAPPED_READ_COUNT_FLAG',
        'MEAN_COVERAGE_FLAG',
        'MEAN_MQ_FLAG'
    ]

    # Convert flag values to scores
    for col in flag_columns:
        df[col + '_SCORE'] = df[col].map(flag_score_map)

    # Total score for each cluster
    df['TOTAL_CONF_SCORE'] = df[[col + '_SCORE' for col in flag_columns]].sum(axis=1)

    # Optionally normalize: score out of 12 (6 flags × max score of 2)
    df['NORMALISED_CONF_SCORE'] = df['TOTAL_CONF_SCORE'] / (2 * len(flag_columns))  # Result: 0 to 1 scale

     #df['qseq_pc_mapping_read'] = df['qseq_pc_mapping_read'].round(1)
    df['qseq_pc_mapping_read'] = df['qseq_pc_mapping_read'].apply(lambda x: float("{:.5f}".format(x)))
    df['qseq_mean_depth'] = df['qseq_mean_depth'].apply(lambda x: float("{:.1f}".format(x)))
    df['NORMALISED_CONF_SCORE'] = df['NORMALISED_CONF_SCORE'].apply(lambda x: float("{:.3f}".format(x)))
    return df


def save_summary(df, sample_name):
    df = df.sort_values("normalised_conf_score", ascending=False)
    output_file = f"{sample_name}_contigs_with_cov_stats.txt"
    df.to_csv(output_file, index=False, sep="\t")
    print(f"Saved final summary to {output_file}")

def parse_mapping_file(mapping_path):
    """Parse tab-delimited file of read_id and reference."""
    mapping = collections.defaultdict(list)
    with open(mapping_path, 'r') as f:
        for line in f:
            read_id, reference = line.strip().split()
            if reference not in mapping[read_id]:  # Manual deduplication
                mapping[read_id].append(reference)
    return mapping

def get_read_lengths(fasta_path):
    """Return dictionary of read_id: length from fasta file."""
    read_lengths = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        read_lengths[record.id] = len(record.seq)
    return read_lengths

def get_reference_lengths(fasta_path):
    """Return dictionary of reference_id: length from consensus fasta."""
    ref_lengths = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        ref_lengths[record.id] = len(record.seq)
    return ref_lengths

def group_lengths_by_reference(mapping, read_lengths):
    """Return reference: list of lengths."""
    reference_lengths = collections.defaultdict(list)
    for read_id, refs in mapping.items():
        length = read_lengths.get(read_id)
        if length:
            # Only use the first reference assigned (or apply a priority rule)
            ref = refs[0] if isinstance(refs, list) else refs
            reference_lengths[ref].append(length)
    return reference_lengths

def plot_coverage_bar(results_dict, output_path="ref_coverage_bar.png"):
    """
    Bar plot showing % of reads ≥80% reference length for each reference.
    """
    refs = list(results_dict.keys())
    fractions = [v["fraction"] for v in results_dict.values()]

    colors = ['green' if v["passes"] else 'red' for v in results_dict.values()]

    plt.figure(figsize=(12, 6))
    bars = plt.bar(refs, fractions, color=colors)
    plt.axhline(y=0.05, color='blue', linestyle='--', label="5% threshold")

    plt.xticks(rotation=90)
    plt.ylim(0, 1)
    plt.ylabel("Fraction of reads ≥ 80% of reference")
    plt.title("Read coverage across references (5/80 rule)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"Bar chart saved as {output_path}")

def analyze_read_lengths_against_reference(reference_lengths, grouped_read_lengths, crl, rpc):
    """
    Analyze read lengths to determine if >=rpc% of reads are >=crl% of the reference length.
    For each reference, check if 5% of reads are ≥80% of reference length
    Parameters:
        read_lengths (dict): read_id -> length of read
        reference_lengths (dict): reference name -> reference length
        crl (int): cutoff percent of reference length
        rpc (int): required percent of reads passing

    Returns:
        df (DataFrame): DataFrame with reference and pass/fail status
    """
    
    results = {}
    for ref, lengths in grouped_read_lengths.items():
        if ref not in reference_lengths:
            print(f"Warning: {ref} not found in consensus fasta.")
            continue

        num_reads = len(lengths)
        ref_len = reference_lengths[ref]
        threshold = (crl / 100) * ref_len

        num_passing = sum(1 for l in lengths if l >= threshold)
        fraction = num_passing / num_reads if num_reads > 0 else 0
        passes = fraction >= (rpc / 100)


        results[ref] = {
            "ref_len": ref_len,
            "num_reads": num_reads,
            "num_passing": num_passing,
            "fraction": fraction,
            "passes": passes
        }

        print(f"{ref}: RefLen={ref_len}, Reads={num_reads}, >={crl}%Ref={num_passing}")
        print(f"Passes {rpc}/{crl}? {'YES' if passes else 'NO'}\n")
    # Save results to file
    df = pd.DataFrame([
    { "qseqid": ref, f"num_passing_{crl}": res["num_passing"], f"pc_read_length_passes_{crl}_{rpc}": res["passes"] }
    for ref, res in results.items()
    ])
    print(df)
    return df

def main():
    args = parse_args()
    blast_df = pd.read_csv(args.blastn_results, sep="\t", header=0)

    filtered_read_counts = parse_bbsplit_log(args.bbsplit_stats)
    print(f"Filtered read counts: {filtered_read_counts}")
    samtools_cov, mosdepth_df, mq_df = load_and_prepare_data(
        args.coverage,
        args.bed,
        args.mapping_quality,
        filtered_read_counts
    )

    
    merged_df = merge_dataframes(samtools_cov, mosdepth_df, mq_df)
    print(merged_df.head())
    #columns_to_extract = ["sacc", "species_updated", "full_lineage" ]
    #subset_df = blast_df[columns_to_extract].drop_duplicates()
    #clean the qseqid by removing version numbers
    #merged_df["qseqid_clean"] = merged_df["qseqid"].str.split(".").str[0]
    #clean the sacc by removing version numbers
    #blast_df["sacc"] = blast_df["sacc"].astype(str).str.strip()
    print(blast_df.head())
    merged_df2 = pd.merge(blast_df, merged_df, left_on="qseqid", right_on="qseqid", how="inner")

    flagged_df = apply_qc_flags(merged_df2)
    flagged_df = flagged_df.rename(columns={
        "species_updated": "taxon_name",
        "query_match_length": "reference_length",
        "qseq_mapping_read_count": "mapping_read_count",
        "qseq_mean_depth": "mean_depth",
        "qseq_pc_mapping_read": "pc_mapping_reads",
        "qseq_pc_cov_30X": "pc_cov_30X",
        "mean_MQ": "mean_mapping_quality",
        "30X_COVERAGE_FLAG": "30x_cov_flag",
        "MAPPED_READ_COUNT_FLAG": "read_count_flag",
        "MEAN_COVERAGE_FLAG": "mean_depth_flag",
        "MEAN_MQ_FLAG": "mean_mq_flag",
        "TOTAL_CONF_SCORE": "total_conf_score",
        "NORMALISED_CONF_SCORE": "normalised_conf_score"
    })
    
    
    flagged_df_subset = flagged_df[["sample_name", "qseqid", "contig_seq", "sacc", "alignment_length", "evalue", "bitscore", "pident", "mismatch",
                     "gapopen", "qstart", "qend", "qlen", "sstart", "send", "slen", "sstrand",
                     "qcovhsp", "staxids", "qseq", "sseq", "qcovs", 
                     "species", "taxon_name", "RNA_type", "stitle", "full_lineage", "ncontigs_per_sacc", "ncontigs_per_spp", "assembly_kmer_cov", "total_score_spp", "total_score_spp_rna", "total_score_sacc", "term_filter", "cov_filter", "best_contig_per_sp_filter", "best_contig_per_sp_rna_filter", "best_contig_per_acc_filter", "mapping_read_count",
                     'pc_mapping_reads', 'mean_depth', 'pc_cov_30X',  'mean_mapping_quality', 'read_count_flag', 'mean_depth_flag', '30x_cov_flag', 'mean_mq_flag', 'total_conf_score','normalised_conf_score']]
 
    
    save_summary(flagged_df_subset, args.sample)

if __name__ == "__main__":
    main()



#ingularity exec -B /work/daff_viral_rnaseq/ ~/.nextflow/NXF_SINGULARITY_CACHEDIR/docker.io-gauthiem-python312.img /work/daff_viral_rnaseq/daff-viral-assembly/bin/derive_coverage_stats.py --sample Ta02_sub --blastn_results /work/daff_viral_rnaseq/daff-viral-assembly/results/Ta02_sub/07_annotation/Ta02_sub_megablast_top_viral_hits_filtered_with_contigs.txt --fastp /work/daff_viral_rnaseq/daff-viral-assembly/results/Ta02_sub/03_fastqc_trimmed/Ta02_sub.fastp.json --coverage /work/daff_viral_rnaseq/daff-viral-assembly/results/Ta02_sub/08_mapping_to_ref/Ta02_sub_coverage.txt  --bed /work/daff_viral_rnaseq/daff-viral-assembly/Ta02_sub.thresholds.bed --consensus /work/daff_viral_rnaseq/daff-viral-assembly/results/Ta02_sub/08_mapping_to_ref/Ta02_sub_bcftools_masked_consensus.fasta --mapping_quality /work/daff_viral_rnaseq/daff-viral-assembly/results/Ta02_sub/08_mapping_to_ref/Ta02_sub_mapq.txt  --reference /work/daff_viral_rnaseq/daff-viral-assembly/results/Ta02_sub/08_mapping_to_ref/Ta02_sub_ref_sequences.fasta