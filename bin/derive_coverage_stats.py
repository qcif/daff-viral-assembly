#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
from functools import reduce
import re


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate coverage statistics for contigs or references"
    )
    parser.add_argument("--mode", choices=["contig", "reference"], required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--blastn_results", required=True)
    parser.add_argument("--bbsplit_stats", required=True)
    parser.add_argument("--coverage", required=True)
    parser.add_argument("--bed", required=True)
    parser.add_argument("--mapping_quality", required=True)
    parser.add_argument("--per_base", required=True)

    return parser.parse_args()

# -----------------------------
# Parse bbsplit stats
# -----------------------------
def parse_bbsplit_log(logfile):

    total_reads = None
    r1_mapped = r2_mapped = 0
    section = None

    with open(logfile) as f:
        for line in f:
            line = line.strip()

            m_total = re.match(r"Reads Used:\s+(\d+)", line)
            if m_total:
                total_reads = int(m_total.group(1))
                continue

            if line.startswith("Read 1 data"):
                section = "R1"
                continue
            if line.startswith("Read 2 data"):
                section = "R2"
                continue

            m_mapped = re.match(r"mapped:\s+[\d\.%]+\s+(\d+)", line)
            if m_mapped:
                count = int(m_mapped.group(1))
                if section == "R1":
                    r1_mapped = count
                elif section == "R2":
                    r2_mapped = count

    if total_reads is None:
        raise ValueError("Could not find total reads in log file")

    total_mapped_reads = r1_mapped + r2_mapped
    clean_reads = total_reads - total_mapped_reads

    return clean_reads

# -----------------------------
# Load coverage stats
# -----------------------------
def load_coverage_stats(coverage_path, bed_path, mq_path, filtered_reads):

    samtools_cov = pd.read_csv(
        coverage_path,
        sep="\t",
        usecols=["#rname", "endpos", "numreads", "meandepth"]
    )

    samtools_cov.rename(
        columns={
            "#rname": "qseqid",
            "endpos": "reference_length",
            "numreads": "mapping_read_count",
            "meandepth": "mean_depth"
        },
        inplace=True
    )
    samtools_cov["pc_mapping_reads"] = (
        samtools_cov["mapping_read_count"] / filtered_reads * 100
    )

    #
    mosdepth = pd.read_csv(bed_path, sep="\t", header=0)
    mosdepth.columns = [
        "qseqid",
        "start",
        "end",
        "region",
        "bases_1x",
        "bases_30x"
    ]

    mosdepth["pc_cov_1X"] = (
                mosdepth["bases_1x"] / mosdepth["end"] * 100
    ).clip(upper=100).round(1)
    mosdepth["pc_cov_30X"] = (
        mosdepth["bases_30x"] / mosdepth["end"] * 100
    ).clip(upper=100).round(1)

    mosdepth = mosdepth[["qseqid", "pc_cov_1X", "pc_cov_30X"]]

    mq = pd.read_csv(mq_path, sep="\t", header=None)
    mq.columns = ["qseqid", "mean_mapping_quality"]

    return samtools_cov, mosdepth, mq

# -----------------------------
# Merge coverage tables
# -----------------------------
def merge_tables(*dfs):

    return reduce(
        lambda left, right: pd.merge(
            left,
            right,
            on="qseqid",
            how="outer"
        ),
        dfs
    ).fillna(0)


# -----------------------------
# QC flags
# -----------------------------
def _apply_shared_qc_flags(df):
    df["30x_cov_flag"] = np.select(
        [
            df.pc_cov_30X >= 100,
            df.pc_cov_30X >= 50,
            df.pc_cov_30X < 50
        ],
        ["GREEN", "ORANGE", "RED"],
        default="GREY"
    )

    df["mean_depth_flag"] = np.select(
        [
            df.mean_depth >= 100,
            df.mean_depth >= 30,
            df.mean_depth < 30
        ],
        ["GREEN", "ORANGE", "RED"],
        default="GREY"
    )

    # df["mean_mq_flag"] = np.select(
    #     [
    #         df.mean_mapping_quality >= 30,
    #         df.mean_mapping_quality >= 10,
    #         df.mean_mapping_quality < 10
    #     ],
    #     ["GREEN", "ORANGE", "RED"],
    #     default="GREY"
    # )

    return df


def _add_confidence_scores(df, flags):
    flag_map = {"GREEN": 2, "ORANGE": 1, "RED": 0, "GREY": 0}

    for f in flags:
        df[f + "_score"] = df[f].map(flag_map)

    df["total_conf_score"] = df[[f + "_score" for f in flags]].sum(axis=1)

    df["normalised_conf_score"] = df["total_conf_score"] / (2 * len(flags))


    return df


def apply_qc_flags_contig(df):
    df = _apply_shared_qc_flags(df)
    alignment_length = pd.to_numeric(df["alignment_length"], errors="coerce")
    subject_length = pd.to_numeric(df["slen"], errors="coerce")
    df["reference_coverage"] = (
        alignment_length.div(subject_length.where(subject_length != 0))
        .mul(100)
        .fillna(0)
    )
    df["reference_coverage_flag"] = np.select(
        [
            df.reference_coverage >= 90,
            df.reference_coverage >= 50,
            df.reference_coverage < 50,
        ],
        ["GREEN", "ORANGE", "RED"],
        default="GREY",
    )

    df["contig_coverage_flag"] = np.select(
        [
            df.qcovs >= 90,
            df.qcovs >= 75,
            df.qcovs < 75,
        ],
        ["GREEN", "ORANGE", "RED"],
        default="GREY",
    )
    
    return _add_confidence_scores(df, [
        "30x_cov_flag",
        "reference_coverage_flag",
        "contig_coverage_flag",
        "mean_depth_flag",
    ])


def apply_qc_flags_reference(df):
    df = _apply_shared_qc_flags(df)
    # df["read_count_flag"] = np.select(
    #     [
    #         df.mapping_read_count >= 1000,
    #         df.mapping_read_count >= 200,
    #         df.mapping_read_count < 200,
    #     ],
    #     ["GREEN", "ORANGE", "RED"],
    #     default="GREY",
    # )
    df["1x_cov_flag"] = np.select(
        [
            df.pc_cov_1X >= 90,
            df.pc_cov_1X >= 80,
            df.pc_cov_1X < 80,
        ],
        ["GREEN", "ORANGE", "RED"],
        default="GREY",
    )
    longest_region = pd.to_numeric(df["longest_covered_region"], errors="coerce")
    reference_length = pd.to_numeric(df["reference_length"], errors="coerce")
    df["pc_longest_covered_region"] = (
        longest_region.div(reference_length.where(reference_length != 0))
        .mul(100).round(1)
        .fillna(0)
    )
    df["1x_long_cov_flag"] = np.select(
        [
            df.pc_longest_covered_region >= 90,
            df.pc_cov_1X >= 40,
            df.pc_cov_1X < 40,
        ],
        ["GREEN", "ORANGE", "RED"],
        default="GREY",
    )
    return _add_confidence_scores(df, [
        "30x_cov_flag",
        "1x_cov_flag",
        "1x_long_cov_flag",
        "mean_depth_flag"
        #"mean_mq_flag",
    ])

# -----------------------------
# Contig mode
# -----------------------------
def process_contigs(blast_df, coverage_df):

    merged = pd.merge(
        blast_df,
        coverage_df,
        on="qseqid",
        how="left"
    ).fillna(0)

    return merged


def longest_covered_regions(depth_file, min_depth=1):
    """Return the longest covered interval length for each reference."""
    longest = {}
    current_ref = None
    run_start = None
    run_end = None

    def save_run():
        if run_start is not None:
            longest[current_ref] = max(
                longest.get(current_ref, 0),
                run_end - run_start,
            )

    with open(depth_file, "r") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue

            ref, start, end, depth = line.split()[:4]
            start, end, depth = int(start), int(end), float(depth)

            if ref != current_ref:
                save_run()
                current_ref = ref
                run_start = run_end = None

            if depth < min_depth:
                save_run()
                run_start = run_end = None
                continue

            if run_start is not None and start <= run_end:
                # Extend across adjacent covered intervals,
                # even when their depths differ.
                run_end = max(run_end, end)
            else:
                save_run()
                run_start, run_end = start, end

        save_run()

    return longest


def longest_covered_regions_dataframe(depth_file, min_depth=1):
    regions = longest_covered_regions(depth_file, min_depth)
    return pd.DataFrame(
        regions.items(),
        columns=["qseqid", "longest_covered_region"],
    )



# -----------------------------
# Reference mode
# -----------------------------
def process_references(blast_df, coverage_df):

    coverage_df["qseqid_clean"] = coverage_df.qseqid.str.split(".").str[0]

    subset = blast_df[
        ["sacc", "species", "species_updated", "full_lineage"]
    ].drop_duplicates()

    subset["sacc"] = subset["sacc"].astype(str).str.strip()

    merged = pd.merge(
        coverage_df,
        subset,
        left_on="qseqid_clean",
        right_on="sacc",
        how="inner"
    )

    return merged


    # -----------------------------
# Save output
# -----------------------------
def save_output(df, sample, mode):
    # Format numeric columns
    df['pc_mapping_reads'] = df['pc_mapping_reads'].round(1)
    df['mean_depth'] = df['mean_depth'].round(1)
    df['normalised_conf_score'] = df['normalised_conf_score'].round(3)
    df = df.sort_values("normalised_conf_score", ascending=False)

    if mode == "contig":
        outfile = f"{sample}_contigs_with_cov_stats.txt"
    else:
        outfile = f"{sample}_reference_with_cov_stats.txt"

    df.to_csv(outfile, sep="\t", index=False)

    print(f"Saved {outfile}")



# -----------------------------
# Main
# -----------------------------
def main():

    args = parse_args()

    blast_df = pd.read_csv(args.blastn_results, sep="\t")

    filtered_reads = parse_bbsplit_log(args.bbsplit_stats)

    samtools_cov, mosdepth_df, mq_df = load_coverage_stats(
        args.coverage,
        args.bed,
        args.mapping_quality,
        filtered_reads
    )

    coverage_df = merge_tables(
        samtools_cov,
        mosdepth_df,
        mq_df,
        longest_covered_regions_dataframe(args.per_base),
    )

    if args.mode == "contig":
        merged = process_contigs(blast_df, coverage_df)
        flagged_df = apply_qc_flags_contig(merged)
        #might want to harmonise between 2 modes here, revisit later
        flagged_df.rename(
            columns={
                "species_updated": "taxon_name"
            },
            inplace=True
        )
        flagged_df_subset = flagged_df[[
            "sample_name","qseqid","contig_seq","sacc","alignment_length","evalue","bitscore","pident","mismatch",
            "gapopen","qstart","qend","qlen","sstart","send","slen","sstrand",
            "qcovhsp","staxids","qseq","sseq","qcovs",
            "species","taxon_name","RNA_type","stitle","full_lineage",
            "ncontigs_per_sacc","ncontigs_per_spp","assembly_kmer_cov",
            "total_score_spp","total_score_spp_rna","total_score_sacc",
            "cov_filter","term_filter",
            "best_contig_per_sp_filter","best_contig_per_sp_rna_filter","best_contig_per_acc_filter",
            "mapping_read_count","pc_mapping_reads","mean_depth","pc_cov_30X",
            "mean_mapping_quality","reference_coverage","reference_coverage_flag","contig_coverage_flag","mean_depth_flag",
            "30x_cov_flag","total_conf_score","normalised_conf_score"
        ]]
        save_output(flagged_df_subset, args.sample, args.mode)

    else:
        merged = process_references(blast_df, coverage_df)
        flagged_df = apply_qc_flags_reference(merged)
        #might want to harmonise between 2 modes here, revisit later
        flagged_df.rename(
            columns={
                "species": "taxon_name",
                "species_updated": "enriched_taxon_name"
            },
            inplace=True
        )
        flagged_df_subset = flagged_df[
            [
                "taxon_name","sacc","full_lineage","reference_length",
                "mapping_read_count","pc_mapping_reads","mean_depth","pc_cov_1X","pc_cov_30X",
                "longest_covered_region","pc_longest_covered_region",
                "mean_mapping_quality","1x_cov_flag","1x_long_cov_flag","mean_depth_flag","30x_cov_flag","total_conf_score","normalised_conf_score"
            ]
        ]

        save_output(flagged_df_subset, args.sample, args.mode)

if __name__ == "__main__":
    main()