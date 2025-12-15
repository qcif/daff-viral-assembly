#!/usr/bin/env python
import argparse
import pandas as pd
import glob
import os
import time
import numpy as np
import re
import json


def read_filtered_read_count(fastp_path):
    with open(fastp_path) as f:
        data = json.load(f)
        raw_reads = data["summary"]["before_filtering"]["total_reads"]
        mean_raw_read_length = data["summary"]["after_filtering"]["read1_mean_length"]
        quality_filtered_reads = data["summary"]["after_filtering"]["total_reads"]
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

    timestr = time.strftime("%Y%m%d-%H%M%S")
    summary_dict = {}

    for raw_read_out in glob.glob("*.fastp.json"):
        mean_raw_read_length = 0
        quality_filtered_reads = 0
        mean_filtered_read_length = 0
        gc_content = 0
        sample = (os.path.basename(raw_read_out).replace('.fastp.json', ''))
        raw_reads, quality_filtered_reads, mean_raw_read_length, mean_filtered_read_length, gc_content  = read_filtered_read_count(raw_read_out)  

        summary_dict[sample] = [raw_reads,  quality_filtered_reads, mean_raw_read_length, mean_filtered_read_length, gc_content]
    print(summary_dict)
        

    for rRNA_clean_read_out in glob.glob("*_rRNA_reads.log"):
        rRNA_cleaned_reads = ()
        sample = (os.path.basename(rRNA_clean_read_out).replace('_rRNA_reads.log', ''))
        with open(rRNA_clean_read_out, 'r') as f:
            for line in f:
                if line.strip().startswith("Result:"):
                    m = re.search(r'Result:\s+(\d+)\s+reads', line)
                    if m:
                        rRNA_cleaned_reads = int(m.group(1))
                        break
            #first_line = next(f)
            #qt_reads = int(first_line[0].strip())
        f.close()
        summary_dict[sample].append(rRNA_cleaned_reads)


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
        phix_cleaned_reads = total_reads - total_mapped_reads
        f.close()
        summary_dict[sample].append(phix_cleaned_reads)


    run_data_df = pd.DataFrame([([k] + v) for k, v in summary_dict.items()], columns=['Sample','raw_reads','quality_filtered_reads','mean_raw_read_length','mean_filtered_read_length','gc_content','rRNA_cleaned_reads','phix_cleaned_reads'])
    run_data_df['percent_qfiltered'] = run_data_df['quality_filtered_reads'] / run_data_df['raw_reads'] * 100
    run_data_df['percent_qfiltered'] = run_data_df['percent_qfiltered'].apply(lambda x: float("{:.2f}".format(x)))
    run_data_df['percent_cleaned'] = run_data_df['phix_cleaned_reads'] / run_data_df['raw_reads'] * 100
    run_data_df['percent_cleaned'] = run_data_df['percent_cleaned'].apply(lambda x: float("{:.2f}".format(x)))
    
    run_data_df = run_data_df.sort_values("Sample")
    run_data_df['raw_reads_flag'] = np.where((run_data_df['raw_reads'] < 8000000), "Less than 8M raw reads", "") # ! confirm with DAFF
    run_data_df['qfiltered_reads_flag'] = np.where((run_data_df['phix_cleaned_reads'] < 2500000), "Less than 2,5M cleaned reads", "") # ! confirm with DAFF
    run_data_df["QC_FLAG"] = np.where(
        (run_data_df['phix_cleaned_reads'] < 2500000), # ! confirm with DAFF
        "RED",
        np.where(
            ((run_data_df['raw_reads'] < 8000000) & # ! confirm with DAFF
            (run_data_df['phix_cleaned_reads'] >= 2500000)), # ! confirm with DAFF
            "ORANGE",
            np.where(
                ((run_data_df['raw_reads'] >= 8000000) & # ! confirm with DAFF 
                (run_data_df['phix_cleaned_reads'] >= 2500000)), # ! confirm with DAFF
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

if __name__ == '__main__':
    main()