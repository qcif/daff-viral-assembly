#!/usr/bin/env python
import sys
import statistics as stats
from collections import defaultdict

pileup_file = sys.argv[1]

# ---- Tunable thresholds ----
SCAN_EDGE = 200
BASELINE_WINDOW = 150
DOWNSTREAM_WINDOW = 20

END_FRAC_THRESH = 0.5
START_FRAC_THRESH = 0.4
DEPTH_DROP_RATIO = 0.6
DEPTH_RISE_RATIO = 1.8
MIN_BASELINE = 20
# ----------------------------

data = defaultdict(list)

# Parse pileup
with open(pileup_file, "r", encoding="latin-1") as f:
    for line in f:
        ctg, pos, ref, dp, bases, qual = line.strip().split("\t")
        pos = int(pos)
        dp = int(dp)

        end_frac = bases.count('$') / dp if dp > 0 else 0
        start_frac = bases.count('^') / dp if dp > 0 else 0

        data[ctg].append((pos, dp, end_frac, start_frac))
        
# print("\n--- SANITY CHECK ---")

# for ctg, rows in data.items():
#     print(f"\nContig: {ctg}")
#     print(f"Number of positions: {len(rows)}")

#     # print first 5
#     print("First 5 entries:")
#     for r in rows[:110]:
#         print(r)

#     # print last 5
#     print("Last 50 entries:")
#     for r in rows[-100:]:
#         print(r)

#     break  # remove this if you want all contigs



for ctg, rows in data.items():

    # Ensure sorted
    rows.sort(key=lambda x: x[0])

    positions = [r[0] for r in rows]
    dp_vals = [r[1] for r in rows]
    end_fracs = [r[2] for r in rows]
    start_fracs = [r[3] for r in rows]

    n = len(rows)
    #print(f"\nProcessing contig {ctg} with {n} positions")

    left_trim = positions[0]
    right_trim = positions[-1]

    # -------------------------
    # 3′ RIGHT END (last 200 bp)
    # -------------------------
    start_scan = max(0, n - SCAN_EDGE)
    for i in range(n-3, start_scan, -1):

        if dp_vals[i] < MIN_BASELINE:
            continue   # ignore zero/very low coverage region

        baseline_start = max(0, i-BASELINE_WINDOW)
        baseline_slice = dp_vals[baseline_start:i]

        if len(baseline_slice) < 20:
            continue

        baseline = stats.median(baseline_slice)

        if baseline < MIN_BASELINE:
            continue

        depth_ratio = dp_vals[i+1] / dp_vals[i] if dp_vals[i] > 0 else 1
        downstream = stats.median(dp_vals[i+1:min(n, i+1+DOWNSTREAM_WINDOW)])

        if (
            end_fracs[i] >= END_FRAC_THRESH and
            depth_ratio < DEPTH_DROP_RATIO and
            downstream < 0.5 * baseline
        ):
            right_trim = positions[i]
            break
    

    # -------------------------
    # 5′ LEFT END (first 200 bp)
    # -------------------------
    end_scan = min(n, SCAN_EDGE)

    for i in range(2, end_scan):

        # Ignore zero / extremely low coverage
        if dp_vals[i] < MIN_BASELINE:
            continue

        # Avoid division by zero
        if dp_vals[i-1] == 0:
            continue

        upstream = stats.median(dp_vals[max(0, i-30):i])
        interior = stats.median(dp_vals[i:min(n, i+BASELINE_WINDOW)])

        depth_ratio = dp_vals[i] / dp_vals[i-1]

        if (
            start_fracs[i] >= START_FRAC_THRESH and
            depth_ratio > DEPTH_RISE_RATIO and
            interior > 1.5 * upstream
        ):
            left_trim = positions[i]
            break

    # Safety guard
    if left_trim >= right_trim:
        left_trim = positions[0]
        right_trim = positions[-1]

    print(ctg, left_trim, right_trim)