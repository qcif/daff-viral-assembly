#!/usr/bin/env python
import argparse
import pandas as pd
import glob
import os
import time
import numpy as np
import re
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def read_filtered_read_count(fastp_path):
    with open(fastp_path) as f:
        data = json.load(f)
        raw_reads = (data["summary"]["before_filtering"]["total_reads"]) // 2  # divide by 2 for paired-end
        mean_raw_read_length = data["summary"]["before_filtering"]["read1_mean_length"]
        quality_filtered_reads = (data["summary"]["after_filtering"]["total_reads"]) // 2  # divide by 2 for paired-end
        mean_filtered_read_length = data["summary"]["after_filtering"]["read1_mean_length"]
        gc_content = data["summary"]["after_filtering"]["gc_content"]
    f.close()
            #raw_read1_section_found = False
        #filtered_read1_section_found = False
        #for line in f:
        #    if "Read1 before filtering" in line:
        #        raw_read1_section_found = True
        #    elif raw_read1_section_found and "total reads:" in line:
        #        # Split on ':' and convert the second part to int
        #        raw_reads = int(line.strip().split(":")[1].strip())
        #    if "Read1 after filtering:" in line:
        #        filtered_read1_section_found = True
        #    elif filtered_read1_section_found and "total reads:" in line:
                # Split on ':' and convert the second part to int
        #        filtered_reads = int(line.strip().split(":")[1].strip())
        #print(f"Filtered reads: {filtered_reads}")
    return raw_reads, quality_filtered_reads,  mean_raw_read_length,  mean_filtered_read_length, gc_content

def main():

    parser = argparse.ArgumentParser(description="Generate QC report for sequencing runs")
    parser.add_argument("--qfiltered_reads_threshold", type=int, default=2500000, help="Threshold for quality filtered reads")
    parser.add_argument("--raw_reads_threshold", type=int, default=8000000, help="Threshold for raw reads")
    args = parser.parse_args()

    timestr = time.strftime("%Y%m%d-%H%M%S")
    summary_dict = {}
    

    for raw_read_out in glob.glob("*.fastp.json"):
        sample = (os.path.basename(raw_read_out).replace('.fastp.json', ''))
        raw_reads, quality_filtered_reads, mean_raw_read_length, mean_filtered_read_length, gc_content  = read_filtered_read_count(raw_read_out)  
        summary_dict[sample] = {
            "raw_reads": raw_reads,
            "quality_filtered_reads": quality_filtered_reads,
            "mean_raw_read_length": mean_raw_read_length,
            "mean_filtered_read_length": mean_filtered_read_length,
            "gc_content": gc_content,
            "rRNA_cleaned_reads": np.nan,
            "phix_cleaned_reads": np.nan,
        }
        #summary_dict[sample] = [raw_reads,  quality_filtered_reads, mean_raw_read_length, mean_filtered_read_length, gc_content]
    print(summary_dict)
        

    for rRNA_clean_read_out in glob.glob("*_bbduk.log"):
        rRNA_cleaned_reads = 0
        sample = (os.path.basename(rRNA_clean_read_out).replace('_bbduk.log', ''))
        with open(rRNA_clean_read_out, 'r') as f:
            for line in f:
                if line.strip().startswith("Result:"):
                    m = re.search(r'Result:\s+(\d+)\s+reads', line)
                    if m:
                        rRNA_cleaned_reads = (int(float(m.group(1)))) // 2  # divide by 2 for paired-end
                        break
            #first_line = next(f)
            #qt_reads = int(first_line[0].strip())
        f.close()
        summary_dict[sample]["rRNA_cleaned_reads"] = rRNA_cleaned_reads


    for bbsplit_log in glob.glob("*_bbsplit_stats.txt"):
        phix_cleaned_reads = ()
        sample = (os.path.basename(bbsplit_log).replace('_bbsplit_stats.txt', ''))
        r1_mapped = r2_mapped = 0

        with open(bbsplit_log, "r") as f:
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
        phix_cleaned_reads = (total_reads - total_mapped_reads) // 2  # divide by 2 for paired-end
        f.close()
        summary_dict[sample]["phix_cleaned_reads"] = phix_cleaned_reads

    run_data_df = (
        pd.DataFrame.from_dict(summary_dict, orient="index")
        .reset_index()
        .rename(columns={"index": "Sample"})
    )
   #run_data_df = pd.DataFrame([([k] + v) for k, v in summary_dict.items()], columns=['Sample','raw_reads','quality_filtered_reads','mean_raw_read_length','mean_filtered_read_length','gc_content','rRNA_cleaned_reads','phix_cleaned_reads'])
    run_data_df['percent_qfiltered'] = run_data_df['quality_filtered_reads'] / run_data_df['raw_reads'] * 100
    run_data_df['percent_qfiltered'] = run_data_df['percent_qfiltered'].apply(lambda x: float("{:.2f}".format(x)))
    run_data_df['percent_cleaned'] = run_data_df['phix_cleaned_reads'] / run_data_df['raw_reads'] * 100
    run_data_df['percent_cleaned'] = run_data_df['percent_cleaned'].apply(lambda x: float("{:.2f}".format(x)))
    int_cols = [
        'raw_reads',
        'quality_filtered_reads',
        'rRNA_cleaned_reads',
        'phix_cleaned_reads'
    ]

    run_data_df[int_cols] = (
        run_data_df[int_cols]
        .apply(pd.to_numeric, errors='coerce')
        .round(0)                 # <-- KEY FIX
        .astype('Int64')
    )
    
    run_data_df = run_data_df.sort_values("Sample")
    run_data_df['raw_reads_flag'] = np.where((run_data_df['raw_reads'] < args.raw_reads_threshold), f"Less than {args.raw_reads_threshold} raw reads", "") # ! confirm with DAFF
    run_data_df['qfiltered_reads_flag'] = np.where((run_data_df['phix_cleaned_reads'] < args.qfiltered_reads_threshold), f"Less than {args.qfiltered_reads_threshold} cleaned reads", "") # ! confirm with DAFF
    run_data_df["QC_FLAG"] = np.where(
        (run_data_df['phix_cleaned_reads'] < args.qfiltered_reads_threshold),
        "RED",
        np.where(
            ((run_data_df['raw_reads'] < args.raw_reads_threshold) &
            (run_data_df['phix_cleaned_reads'] >= args.qfiltered_reads_threshold)),
            "ORANGE",
            np.where(
                ((run_data_df['raw_reads'] >= args.raw_reads_threshold) &
                (run_data_df['phix_cleaned_reads'] >= args.qfiltered_reads_threshold)), 
                "GREEN",
                ""
            )
        )
    )

    run_data_df.to_csv("run_qc_report_" + timestr + ".txt", index = None, sep="\t")

    summary_table = run_data_df.to_html(index=False).replace('<table border="1" class="dataframe">','<table class="table table-striped">') # use bootstrap styling
    html_string = '''
    <html>
        <head>
            <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.1/css/bootstrap.min.css">
            <style>body{ margin:0 100; background:whitesmoke; }</style>
        </head>
        <body>
            <h1>Quality check report</h1>
            <blockquote>
            <p><b> Definitions:</b></p>
            <p><b> raw_reads </b>= total raw reads sequenced </p>
            <p><b> processed_reads </b>= total reads left-over after adapter trimming and quality filtering </p>
            <p><b> percent_qfiltered </b>= (quality filtered/raw_reads x 100)</p>
            <p><b> percent_cleaned </b>= (cleaned reads/raw_reads x 100)</p>
            </blockquote>
            <!-- *** Section 1 *** --->
            ''' + summary_table + '''
        </body>
    </html>'''

    report = open("run_qc_report_" + timestr + ".html", "w")
    report.write(html_string)
    report.close()
    

    # === STACKED BAR PLOT (PERCENTAGES, SAMPLE ON Y-AXIS) ===

    # Convert Int64 (nullable) to plain float to avoid index-alignment NaN issues
    _raw   = run_data_df["raw_reads"].astype(float).values
    _qf    = run_data_df["quality_filtered_reads"].astype(float).values
    _rrna  = pd.to_numeric(run_data_df["rRNA_cleaned_reads"], errors="coerce").values  # NaN if step missing
    _phix  = run_data_df["phix_cleaned_reads"].astype(float).values
    _samples = run_data_df["Sample"].values

    read_origin_df = pd.DataFrame({
        "raw_reads-quality_filtered_reads": _raw - _qf,
        "quality_filtered_reads-rRNA_cleaned_reads": np.where(np.isnan(_rrna), 0.0, _qf - _rrna),
        "rRNA_cleaned_reads-phix_cleaned_reads": np.where(np.isnan(_rrna) | np.isnan(_phix), 0.0, _rrna - _phix),
        "phix_cleaned_reads": _phix,
    }, index=_samples)
    read_origin_df = read_origin_df.fillna(0)


    pc_df = read_origin_df.apply(lambda x: 100 * x / float(x.sum()), axis=1)
    print(pc_df)

    # Plot as horizontal stacked bar
    ax = pc_df.plot.barh(
        stacked=True,
        figsize=(10, 7),
        color=["#4C99E7", "#FFD700", "#F37736", "#7BC043" ],
        width=0.75
    )
    ax.set_xlabel("Read Percentage (%)")
    ax.set_ylabel("Sample")
    ax.set_title("Read origin summary per sample (percentage)")
    vline = ax.axvline(x=50, color="black", linestyle="--", linewidth=1)
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    legend_handles = [
        Patch(facecolor="#4C99E7", label="Quality-filtered reads"),
        Patch(facecolor="#FFD700", label="rRNA-matched reads"),
        Patch(facecolor="#F37736", label="PhiX-matched reads"),
        Patch(facecolor="#7BC043", label="Cleaned reads"),
        Line2D([0], [0], color="black", linestyle="--", linewidth=1, label="50% threshold"),
    ]
    ax.legend(handles=legend_handles, bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(f"read_origin_summary.png", dpi=150, bbox_inches='tight')
    plt.savefig(f"read_origin_summary.pdf", bbox_inches='tight')
    plt.close()

    # === PER-SAMPLE BAR PLOTS (ONE PDF PER SAMPLE) ===

    for sample in read_origin_df.index:

        sample_data = read_origin_df.loc[sample]

        # Convert to percentages (safe)
        total = sample_data.sum()
        if total == 0:
            continue

        sample_pct = (sample_data / total) * 100

        # Plot (vertical bar for clarity)
        fig, ax = plt.subplots(figsize=(6, 1.5))
        ax.legend(
            handles=legend_handles,
            loc='center left',
            bbox_to_anchor=(1.0, 0.5),
            fontsize=6,
            frameon=False,
            handlelength=1,
            labelspacing=0.25,
            borderpad=0.2
        )

        # Smaller fonts everywhere
        #ax.set_title(f"{sample}", fontsize=9)
        #ax.set_xlabel("Read %", fontsize=8)

        ax.tick_params(axis='x', labelsize=7)
        ax.tick_params(axis='y', labelsize=7)

        # Define order + labels (IMPORTANT for consistency)
        labels = [
            "Filtered out",
            "rRNA",
            "PhiX",
            "Cleaned reads"
        ]

        values = [
            sample_pct["raw_reads-quality_filtered_reads"],
            sample_pct["quality_filtered_reads-rRNA_cleaned_reads"],
            sample_pct["rRNA_cleaned_reads-phix_cleaned_reads"],
            sample_pct["phix_cleaned_reads"]
        ]

        colors = ["#4C99E7", "#FFD700", "#F37736", "#7BC043"]

        # Build stacked bar manually
        left = 0
        for val, color in zip(values, colors):
            ax.barh(sample, val, left=left, color=color)
            left += val

        ax.set_xlim(0, 100)
        ax.set_xlabel("Read origin (%)")

        # Add 50% line
        ax.axvline(x=50, color="black", linestyle="--", linewidth=1)

        # Legend
        from matplotlib.patches import Patch
        from matplotlib.lines import Line2D

        legend_handles = [
            Patch(facecolor="#4C99E7", label="Quality-filtered reads"),
            Patch(facecolor="#FFD700", label="rRNA-matched reads"),
            Patch(facecolor="#F37736", label="PhiX-matched reads"),
            Patch(facecolor="#7BC043", label="Cleaned reads"),
            Line2D([0], [0], color="black", linestyle="--", linewidth=1, label="50% threshold"),
        ]
        # Remove box, keep axes
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        ax.spines['left'].set_color('#666666')
        ax.spines['bottom'].set_color('#666666')

        ax.spines['left'].set_linewidth(0.8)
        ax.spines['bottom'].set_linewidth(0.8)

        ax.tick_params(axis='both', labelsize=7, length=3, width=0.8, color='#666666')

        # Clean y-axis (only one label)
        ax.set_ylabel("")

        plt.tight_layout(rect=[0, 0, 0.85, 0.9])

        # Save
        plt.savefig(f"{sample}_read_origin_stacked_pct.png", dpi=150, bbox_inches='tight')
        plt.savefig(f"{sample}_read_origin_stacked_pct.pdf", bbox_inches='tight')

        plt.close()
if __name__ == '__main__':
    main()