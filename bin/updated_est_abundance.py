#!/usr/bin/env python
import os, sys, argparse
from time import gmtime, strftime


class Tree:
    def __init__(self, name, taxid, level_num, level_id, all_reads, lvl_reads, children=None, parent=None):
        self.name = name
        self.taxid = taxid
        self.level_num = level_num
        self.level_id = level_id
        self.all_reads = all_reads
        self.lvl_reads = lvl_reads
        self.children = []
        self.parent = parent
        if children:
            for child in children:
                self.add_child(child)

    def add_child(self, node):
        self.children.append(node)


def process_kraken_line(line):
    parts = line.rstrip().split("\t")
    if len(parts) < 7:
        return None
    percent, reads_all, reads_level, bp, extra, level_id, taxid, *name_parts = parts
    name = " ".join(name_parts).strip()
    level_num = int(extra) if extra.isdigit() else 0
    return [name, taxid, level_num, level_id, int(reads_all), int(reads_level)]


def check_report_file(file):
    with open(file) as f:
        first_line = f.readline()
        if len(first_line.split("\t")) < 7:
            sys.stderr.write("ERROR: Kraken report file not in expected format\n")
            exit(1)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True, help='Kraken report file')
    parser.add_argument('-k', '--kmer_distr', required=True, help='Kmer distribution file')
    parser.add_argument('-o', '--output', required=True, help='Bracken output abundance file')
    parser.add_argument('-l', '--level', default='S', help='Desired classification level')
    parser.add_argument('--out-report', default='', help='Optional new Kraken-style report')
    parser.add_argument('-t', '--thresh', default=10, type=int, help='Minimum reads threshold')
    args = parser.parse_args()

    sys.stdout.write("PROGRAM START TIME: " + strftime("%m-%d-%Y %H:%M:%S", gmtime()) + "\n")
    check_report_file(args.input)

    # ------------------------------
    # Containers
    # ------------------------------
    lvl_taxids = {}
    map2lvl_taxids = {}

    lvl_taxids_S = {}
    lvl_taxids_S1 = {}
    map2lvl_taxids_S = {}
    map2lvl_taxids_S1 = {}

    taxid_to_name = {}
    taxid_parent = {}
    s1_parents = set()

    root_node = None
    prev_node = None
    last_seen_S = None
    s1_parents = set()

    # ------------------------------
    # Parse Kraken report
    # ------------------------------
    with open(args.input) as f:
        for line in f:
            if line.startswith(("#", "%")) or len(line.strip()) == 0:
                continue

            parsed = process_kraken_line(line)
            if not parsed:
                continue

            name, taxid, level_num, level_id, all_reads, lvl_reads = parsed

            # Track names
            if taxid not in taxid_to_name:
                taxid_to_name[taxid] = {}
            taxid_to_name[taxid][level_id] = name

            if "S1" in taxid_to_name[taxid]:
                preferred_name = taxid_to_name[taxid]["S1"]
            elif "S" in taxid_to_name[taxid]:
                preferred_name = taxid_to_name[taxid]["S"]
            else:
                preferred_name = name

            if level_id == "U":
                continue

            if taxid == '1':
                root_node = Tree(preferred_name, taxid, level_num, 'R', all_reads, lvl_reads)
                prev_node = root_node
                continue

            if prev_node is None:
                prev_node = root_node

            while prev_node and level_num != prev_node.level_num + 1:
                prev_node = prev_node.parent
            if prev_node is None:
                prev_node = root_node

            curr_node = Tree(preferred_name, taxid, level_num, level_id,
                             all_reads, lvl_reads, parent=prev_node)
            prev_node.add_child(curr_node)
            if level_id == "S":
                last_seen_S = taxid
            elif level_id == "S1" and last_seen_S:
                s1_parents.add(last_seen_S)
            else:
                last_seen_S = None

            prev_node = curr_node

            # ------------------------------
            # Capture S and S1 separately
            # ------------------------------
            if args.level == "S" and all_reads >= args.thresh:
                if level_id == "S":
                    lvl_taxids_S[taxid] = [preferred_name, all_reads, lvl_reads, 0]
                    map2lvl_taxids_S[taxid] = [taxid, lvl_reads, 0]
                elif level_id == "S1":
                    lvl_taxids_S1[taxid] = [preferred_name, all_reads, lvl_reads, 0]
                    map2lvl_taxids_S1[taxid] = [taxid, lvl_reads, 0]

            elif level_id == args.level and all_reads >= args.thresh:
                lvl_taxids[taxid] = [preferred_name, all_reads, lvl_reads, 0]
                map2lvl_taxids[taxid] = [taxid, lvl_reads, 0]
        print("Suppressed S taxids:", s1_parents, file=sys.stderr)

    # ------------------------------
    # Final merge for S/S1 logic
    # ------------------------------
    if args.level == "S":
        lvl_taxids = dict(lvl_taxids_S1)
        map2lvl_taxids = dict(map2lvl_taxids_S1)

        for taxid in lvl_taxids_S:
            if taxid not in s1_parents:
                lvl_taxids[taxid] = lvl_taxids_S[taxid]
                map2lvl_taxids[taxid] = map2lvl_taxids_S[taxid]

    # ------------------------------
    # Read kmer distribution
    # ------------------------------
    kmer_distr_dict = {}
    with open(args.kmer_distr) as kf:
        for line in kf.readlines()[1:]:
            parts = line.strip().split("\t")
            mapped_taxid = parts[0]
            temp_dict = {}
            for g in parts[1].split(" "):
                g_taxid, mkmers, tkmers = g.split(":")
                frac = float(mkmers) / float(tkmers)
                temp_dict[g_taxid] = [frac]
            kmer_distr_dict[mapped_taxid] = temp_dict

    # ------------------------------
    # Distribute reads from parents
    # ------------------------------
    curr_nodes = [root_node]
    while curr_nodes:
        node = curr_nodes.pop(0)
        for child in node.children:
            curr_nodes.append(child)

        if node.level_id == args.level or (node.level_id == "S1" and args.level == "S"):
            continue
        if node.lvl_reads == 0 or node.taxid not in kmer_distr_dict:
            continue

        curr_dict = kmer_distr_dict[node.taxid]
        all_genome_reads = sum(map2lvl_taxids[g][1] if g in map2lvl_taxids else 0 for g in curr_dict)
        if all_genome_reads == 0:
            continue

        for genome in curr_dict:
            if genome not in map2lvl_taxids:
                continue
            add_reads = (curr_dict[genome][0] *
                         map2lvl_taxids[genome][1] /
                         all_genome_reads * node.lvl_reads)
            map2lvl_taxids[genome][2] += add_reads

    # ------------------------------
    # Summarise reads
    # ------------------------------
    for genome in map2lvl_taxids:
        lvl_taxids[genome][3] += map2lvl_taxids[genome][2]

    sum_all_reads = sum(float(all_reads + added)
                        for _, all_reads, _, added in lvl_taxids.values())

    # ------------------------------
    # Write output
    # ------------------------------
    with open(args.output, 'w') as out_f:
        out_f.write('name\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\tadded_reads\tnew_est_reads\tfraction_total_reads\n')
        for taxid, (name, all_reads, lvl_reads, added) in lvl_taxids.items():
            new_all = float(all_reads) + float(added)
            out_f.write(f"{name}\t{taxid}\t{args.level}\t{int(all_reads)}\t"
                        f"{int(new_all - int(all_reads))}\t{int(new_all)}\t"
                        f"{new_all / sum_all_reads:.5f}\n")

    sys.stdout.write(f"BRACKEN OUTPUT PRODUCED: {args.output}\n")
    sys.stdout.write("PROGRAM END TIME: " + strftime("%m-%d-%Y %H:%M:%S", gmtime()) + "\n")


if __name__ == "__main__":
    main()
