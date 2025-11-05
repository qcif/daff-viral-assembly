#!/usr/bin/env python
import argparse
import pandas as pd
import glob
import os
import time
import numpy as np

def main():

    timestr = time.strftime("%Y%m%d-%H%M%S")
    summary_dict = {}

    for raw_read_out in glob.glob("*raw_NanoStats.txt"):
        raw_reads = ()
        line_number = 0
        sample = (os.path.basename(raw_read_out).replace('_raw_NanoStats.txt', ''))
        with open(raw_read_out, 'r') as f:
            for line in f:
                string_to_search = ("number_of_reads")
                line_number += 1
                if string_to_search in line:
                    elements = line.rsplit('\t')
                    raw_reads = int(elements[1].strip())
                    print(raw_reads)

        summary_dict[sample] = [raw_reads]
        print(summary_dict)
        f.close()

    for qt_read_out in glob.glob("*filtered_NanoStats.txt"):
        qt_reads = ()
        line_number = 0 
        sample = (os.path.basename(qt_read_out).replace('_filtered_NanoStats.txt', ''))
        with open(qt_read_out, 'r') as f:
            for line in f:
                string_to_search = ("number_of_reads")
                line_number += 1
                if string_to_search in line:
                    elements = line.rsplit('\t')
                    qt_reads = int(elements[1].strip())
                    print(qt_reads)
            #first_line = next(f)
            #qt_reads = int(first_line[0].strip())
        f.close()
        summary_dict[sample].append(qt_reads)

    run_data_df = pd.DataFrame([([k] + v) for k, v in summary_dict.items()], columns=['Sample','raw_reads','processed_reads'])
    run_data_df['percent_processed'] = run_data_df['processed_reads'] / run_data_df['raw_reads'] * 100
    run_data_df['percent_processed'] = run_data_df['percent_processed'].apply(lambda x: float("{:.2f}".format(x)))
    run_data_df = run_data_df.sort_values("Sample")
    run_data_df['raw_reads_flag'] = np.where((run_data_df['raw_reads'] < 2500), "Less than 2500 raw reads", "")
    run_data_df['processed_reads_flag'] = np.where((run_data_df['processed_reads'] < 200), "Less than 200 processed reads", "")
    run_data_df["QC_FLAG"] = np.where(
        (run_data_df['processed_reads'] < 200),
        "RED",
        np.where(
            ((run_data_df['raw_reads'] < 2500) &
            (run_data_df['processed_reads'] >= 200)),
            "ORANGE",
            np.where(
                ((run_data_df['raw_reads'] >= 2500) &
                (run_data_df['processed_reads'] >= 200)),
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
            <p><b> processed_reads </b>= total reads left-over after adapter trimming and/or quality filtering </p>
            <p><b> percent_processed </b>= (processed_reads/raw_reads x 100)</p>
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