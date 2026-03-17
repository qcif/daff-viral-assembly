#!/usr/bin/env python

import argparse
from time import gmtime, strftime

MAIN_RANKS = ["R","K","D","P","C","O","F","G","S"]
VIRUS_TAXID = "10239"
VIROID_TAXID = "12884"
VIRUS_FRACTION_MIN = 0.000001

class Tree:
    def __init__(self, name, taxid, level_num, level_id, all_reads, lvl_reads, parent=None):
        self.name = name
        self.taxid = taxid
        self.level_num = level_num
        self.level_id = level_id
        self.all_reads = all_reads
        self.lvl_reads = lvl_reads
        self.children = []
        self.parent = parent

    def add_child(self, node):
        self.children.append(node)
def genus_has_unresolved_reads(genus):

    child_sum = sum(child.all_reads for child in genus.children)

    return genus.all_reads > child_sum

def family_has_unresolved_reads(fam):

    child_sum = sum(child.all_reads for child in fam.children)

    return fam.all_reads > child_sum
    
def prefer_strain_node(node):

    s1_node = None

    for child in node.children:
        if child.level_id.startswith("S1"):
            s1_node = child

    return s1_node


def is_viral_node(node):

    curr = node

    while curr:

        if curr.taxid == VIRUS_TAXID or curr.taxid == VIROID_TAXID:
            return True

        curr = curr.parent

    return False

def find_family_node(node):

    curr = node

    while curr:

        if get_base_rank(curr.level_id) == "F":
            return curr

        curr = curr.parent

    return None

def process_kraken_report(line):

    parts = line.rstrip().split("\t")
    if len(parts) < 6:
        return None

    try:
        all_reads = int(parts[1])
        level_reads = int(parts[2])
    except:
        return None

    level_id = parts[-3]
    taxid = parts[-2]
    name = parts[-1]

    spaces = len(name) - len(name.lstrip())
    level_num = spaces // 2

    name = name.strip()

    return name, taxid, level_num, level_id, all_reads, level_reads


def get_base_rank(level):

    if not level:
        return ""

    # Treat strain branches as species
    if level.startswith("S"):
        return "S"

    if level.startswith("G"):
        return "G"

    if level.startswith("F"):
        return "F"

    return level[0]


def rank_name(level):

    base = get_base_rank(level)

    mapping = {
        "S": "species",
        "G": "genus",
        "F": "family",
        "O": "order",
        "C": "class",
        "P": "phylum",
        "D": "domain",
        "K": "kingdom"
    }

    if base not in mapping:
        return "unclassified"

    # Treat strain levels (S1,S2...) as species
    if base == "S":
        return "species"

    if len(level) > 1:
        return mapping[base] + "_unclassified"

    return mapping[base]


def has_deeper_rank(node, rank_index):

    for child in node.children:

        base = get_base_rank(child.level_id)

        if base in MAIN_RANKS:

            if MAIN_RANKS.index(base) > rank_index:
                return True

        if has_deeper_rank(child, rank_index):
            return True

    return False

def collect_lowest_nodes(node, thresh):

    selected = []

    # Traverse children first
    for child in node.children:
        selected.extend(collect_lowest_nodes(child, thresh))

    # Handle explicit "unclassified" nodes but do NOT block children
    if node.name.lower().startswith("unclassified"):
        if node.all_reads >= thresh:
            selected.append(node)
        return selected

    base = get_base_rank(node.level_id)

    if base not in MAIN_RANKS:
        return selected

    rank_index = MAIN_RANKS.index(base)

    # Prefer S1 strain nodes
    if base == "S":
        s1 = prefer_strain_node(node)
        if s1 and s1.all_reads >= thresh:
            return selected + [s1]

    # Prevent S2 replacing S1
    if node.level_id.startswith("S2"):
        if node.parent and node.parent.level_id.startswith("S1"):
            return selected

    # Standard lowest-node selection
    if node.all_reads >= thresh and not has_deeper_rank(node, rank_index):
        selected.append(node)

    # -------- FAMILY RESCUE --------
    # -------- FAMILY RESCUE --------
    if is_viral_node(node):

        fam = find_family_node(node)

        if fam and fam.all_reads >= thresh:

            if family_has_unresolved_reads(fam) and fam not in selected:
                selected.append(fam)
    # -------- GENUS RESCUE --------
    if is_viral_node(node):

        if get_base_rank(node.level_id) == "G":

            if node.all_reads >= thresh:

                if genus_has_unresolved_reads(node) and node not in selected:
                    selected.append(node)
    return selected

def build_tree(report):

    root = None
    prev_node = None

    for line in report:

        vals = process_kraken_report(line)

        if not vals:
            continue

        name, taxid, level_num, level_id, all_reads, lvl_reads = vals

        # skip unclassified rows
        if level_id == "U":
            continue

        # root node
        if taxid == "1":
            root = Tree(name, taxid, level_num, "R", all_reads, lvl_reads)
            prev_node = root
            continue

        # skip anything before root
        if root is None:
            continue

        # move up until correct parent
        while prev_node and level_num != prev_node.level_num + 1:
            prev_node = prev_node.parent

        # if traversal fell off tree, reset to root
        if prev_node is None:
            prev_node = root

        node = Tree(name, taxid, level_num, level_id, all_reads, lvl_reads, prev_node)

        prev_node.add_child(node)
        prev_node = node

    return root


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("-i","--input", required=True)
    parser.add_argument("-o","--output", required=True)
    parser.add_argument("-t","--thresh", default=10, type=int)

    args = parser.parse_args()

    print("PROGRAM START TIME:", strftime("%Y-%m-%d %H:%M:%S", gmtime()))

    with open(args.input) as f:
        root = build_tree(f)

    targets = collect_lowest_nodes(root, args.thresh)

    lvl_taxids = {}

    total_reads = 0

    for node in targets:
        #Use lvl_reads for explicit "unclassified" nodes, and keep all_reads for normal retained nodes like S1.
        #this cactehes cases where an "unclassified" node is retained instead of deeper taxa, 
        # and prevents double counting of reads when both an "unclassified" node and deeper taxa are retained for the same branch.
        # Avoid double counting when parent and children both retained
        if "unclassified" in node.name.lower():
            reads_to_report = node.lvl_reads

        elif get_base_rank(node.level_id) in ["F", "G"] and node.children:
            reads_to_report = node.lvl_reads

        else:
            reads_to_report = node.all_reads


        base_rank = get_base_rank(node.level_id)
        lowest = rank_name(node.level_id)

        # If family node has unresolved reads, mark as family_unclassified
        #If a family node still has reads assigned directly at that level (lvl_reads > 0) while deeper taxa exist, 
        # report it as family_unclassified.
        if base_rank == "F" and node.lvl_reads > 0 and node.children:
            lowest = "family_unclassified"

        elif base_rank == "G" and node.lvl_reads > 0 and node.children:
            lowest = "genus_unclassified"

        lvl_taxids[node.taxid] = {
        "name": node.name,
        "rank": base_rank,
        "lowest_rank": lowest,
        "kraken_reads": reads_to_report,
}
        total_reads += reads_to_report


    with open(args.output,"w") as out:

        out.write(
            "name\ttaxonomy_id\ttaxonomy_lvl\tlowest_rank\treads\tfraction_total_reads\n"
        )

        for taxid, data in lvl_taxids.items():

            reads = data["kraken_reads"]

            frac = 0
            if total_reads > 0:
                frac = reads / total_reads
            if frac < VIRUS_FRACTION_MIN:
                continue
            out.write(
                f"{data['name']}\t{taxid}\t{data['rank']}\t{data['lowest_rank']}\t"
                f"{reads}\t{frac:.6f}\n"
            )

    print("OUTPUT WRITTEN:", args.output)
    print("PROGRAM END TIME:", strftime("%Y-%m-%d %H:%M:%S", gmtime()))


if __name__ == "__main__":
    main()