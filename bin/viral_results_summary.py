#!/usr/bin/env python
import argparse
import pandas as pd
import re
import os.path
from Bio import SeqIO
from analyses_config import file_exists, is_file_empty

# Define PFAMs of interest
PFAMS = [
    "PF00863.26","PF13608.11","PF00998.29","PF00680.26","PF00078.33","PF00767.23",
    "PF08440.16","PF07184.16","PF05533.19","PF01660.23","PF01443.25","PF17646.7",
    "PF01785.23","PF05520.16","PF06922.16","PF03225.20","PF05733.17","PF11479.14",
    "PF11757.14","PF04808.17","PF19223.6","PF00978.27","PF05379.18","PF01307.23",
    "PF20896.3","PF00721.28","PF01107.24","PF01577.23","PF07652.21","PF12524.14", 
    "PF02123.22", "PF00894.25", "PF00729.25", "PF02122.22"
]

# PFAM renaming dictionary
PFAM_RENAME = {
    "PF00863.26": "PF00863.26 (Peptidase_S8_S53)",
    "PF13608.11": "PF13608.11 (Protein P3 of Potyviral polyprotein)",
    "PF00998.29": "PF00998.29 (Viral RNA dependent RNA polymerase)",
    "PF00680.26": "PF00680.26 (Viral RNA-dependent RNA polymerase)",
    "PF00078.33": "PF00078.33 (Reverse transcriptase (RNA-dependent DNA polymerase))",
    "PF00767.23": "PF00767.23 (Potyvirus coat protein)",
    "PF08440.16": "PF08440.16 (Potyviridae polyprotein)",
    "PF07184.16": "PF07184.16 (Citrus tristeza virus P33 protein)",
    "PF05533.19": "PF05533.19 (Peptidase_C42)",
    "PF01660.23": "PF01660.23 (Viral methyltransferase)",
    "PF01443.25": "PF01443.25 (Viral superfamily 1 RNA helicase core domain)",
    "PF17646.7": "PF17646.7 (Closterovirus 1a polyprotein central region)",
    "PF01785.23": "PF01785.23 (Closterovirus coat protein)",
    "PF05520.16": "PF05520.16 (Citrus_P18)",
    "PF06922.16": "PF06922.16 (Citrus tristeza virus P13 protein)",
    "PF03225.20": "PF03225.20 (Viral heat shock protein Hsp90 homologue)",
    "PF05733.17": "PF05733.17 (Tenuivirus/Phlebovirus nucleocapsid protein)",
    "PF11479.14": "PF11479.14 (RNA silencing suppressor P21 C-terminal domain)",
    "PF11757.14": "PF11757.14 (Suppressor of RNA silencing P21-like N-terminal domain)",
    "PF04808.17": "PF04808.17 (Citrus tristeza virus (CTV) P23 protein)",
    "PF19223.6": "PF19223.6 (Chroparavirus methyltransferase)",
    "PF00978.27": "PF00978.27 (RNA dependent RNA polymerase)",
    "PF05379.18": "PF05379.18 (Carlavirus endopeptidase)",
    "PF01307.23": "PF01307.23 (Plant viral movement protein)",
    "PF20896.3": "PF20896.3 (Tomato mosaic virus helicase, N-terminal domain)",
    "PF00721.28": "PF00721.28 (Virus coat protein (TMV like))",
    "PF01107.24": "PF01107.24 (Viral movement protein (MP))",
    "PF01577.23": "PF01577.23 (Potyvirus P1 protease)",
    "PF07652.21": "PF07652.21 (Flavivirus DEAD domain)",
    "PF12524.14": "PF12524.14 (sDNA virus glycoprotein L C terminal)",
    "PF02123.22": "PF02123.22 (Viral RNA-directed RNA-polymerase)",
    "PF00894.25": "PF00894.25 (Luteovirus coat protein)",
    "PF00729.25": "PF00729.25 (Viral coat protein (S domain))",
    "PF02122.22": "PF02122.22 (Peptidase S39)"
}

RDRP_PFAMS = ["PF00998.29","PF00680.26","PF00078.33","PF00978.27","PF02123.22"]

def parse_hmmscan_per_target(path):
    rows = []
    with open(path, "r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = re.split(r"\s+", line.strip(), maxsplit=18)
            rows.append(parts)

    columns = [
        "target_name","target_accession","query_name","query_accession",
        "full_sequence_evalue","full_sequence_score","full_sequence_bias",
        "best_domain_evalue","best_domain_score","best_domain_bias",
        "domain_n_exp","domain_n_reg","domain_n_clu","domain_n_ov",
        "domain_n_env","domain_n_dom","domain_n_rep","domain_n_inc",
        "description_of_target"
    ]

    df = pd.DataFrame(rows, columns=columns)
    df["full_sequence_score"] = df["full_sequence_score"].astype(float)
    return df

def filter_hmmscan(df):
    df_filtered = df[df["full_sequence_score"] >= 50.0]
    df_filtered = df_filtered[df_filtered["target_accession"].isin(PFAMS)]
    df_filtered = df_filtered.drop_duplicates(subset=["query_name","target_accession"])
    return df_filtered

def summarize_domains(df_filtered):
    # Extract ORF number and NODE name
    df_filtered["ORFs"] = df_filtered["query_name"].str.extract(r"_ORF\.(\d+)$")
    df_filtered["query_name"] = df_filtered["query_name"].str.replace(r"_ORF\.\d+$", "", regex=True)

    # Determine which PFAM accessions are present
    pfams_present = sorted(df_filtered["target_accession"].unique())

    # Create count columns dynamically
    for pf in pfams_present:
        df_filtered[pf] = df_filtered["target_accession"].apply(lambda x: 1 if x == pf else 0)

    # Build aggregation dictionary
    agg_dict = {"ORFs": lambda x: ",".join(sorted(x.dropna(), key=lambda v: int(v)))}
    for pf in pfams_present:
        agg_dict[pf] = "sum"

    # Group by NODE
    summary = df_filtered.groupby("query_name").agg(agg_dict).reset_index()

    # Add Total column (sum of PFAM counts)
    summary["PFAM_total"] = summary[pfams_present].sum(axis=1)

    # Remove rows where PFAM count is zero
    summary = summary[summary["PFAM_total"] > 0]

    # If no domains detected → return empty table (no totals row)
    if summary.empty:
        base_cols = ["query_name", "ORFs", "PFAM_total", "RdRp"]
        return pd.DataFrame(columns=base_cols)
    
    # Add RdRp column
    existing_rdRp_pfams = [pf for pf in RDRP_PFAMS if pf in summary.columns]
    #summary["RdRp"] = summary[existing_rdRp_pfams].sum(axis=1) > 0
    summary["RdRp"] = summary[existing_rdRp_pfams].sum(axis=1) > 0 if existing_rdRp_pfams else False


    # Sort by highest number of PFAM occurrences
    summary = summary.sort_values(by="PFAM_total", ascending=False)

    # Rename PFAMs if mapping exists
    rename_dict = {pf: PFAM_RENAME[pf] for pf in pfams_present if pf in PFAM_RENAME}
    summary.rename(columns=rename_dict, inplace=True)

    # Build total summary column
    pfam_cols_final = [rename_dict.get(pf, pf) for pf in pfams_present]
    total_row = summary[pfam_cols_final].sum()
    total_row["PFAM_total"] = summary["PFAM_total"].sum()
    #fill the column row with adequate information
    #total_row["query_name"] = "Total_counts"
    #total_row["ORFs"] = ""
    total_row["Total_counts"] = ""
    #summary = pd.concat([summary, total_row.to_frame().T], ignore_index=True)

    #add column Total
    summary = pd.concat(
        [summary, total_row.to_frame().T],
        ignore_index=True
    )
    return summary
def add_group_max(df, group_col, target_col, new_col):
    df[new_col] = df.groupby(group_col)[target_col].transform("max")
    return df

VIROID_RENAME = {
    "Citrus dwarfing viroid": "Apscaviroid nanocitri",
    "Citrus viroid IIIa": "Apscaviroid nanocitri"
}

def standardize_viroid_names(df, column, df_name=""):
    """Ensure consistent naming for two viroids across any dataframe.

    Converts both known variants to the canonical name and warns if any
    old names remain afterwards.
    """
    if column not in df.columns:
        return df

    #df[column] = df[column].replace(VIROID_RENAME)
    df = df.copy()
    df.loc[:, column] = df[column].replace(VIROID_RENAME)

    remaining = df[column].isin(VIROID_RENAME.keys())
    if remaining.any():
        remaining_vals = df.loc[remaining, column].unique().tolist()
        print(f"WARNING: {df_name} still contains outdated viroid names in '{column}': {remaining_vals}")

    return df

# Function to convert FASTA file to DataFrame
def fasta_to_dataframe(fasta_file):
    records = SeqIO.parse(fasta_file, "fasta")
    
    # List to hold the sequence data
    data = []
    
    for record in records:
        # Append ID and sequence to the list
        data.append([record.id, str(record.seq)])
    
    # Create a DataFrame
    df = pd.DataFrame(data, columns=["seq_name", "contig_seq"])
    return df


def enrich_with_taxonomy(df, taxonomy):
    """Add family from the taxonomy file using accession (version-agnostic)."""

    if df.empty:
        df["family"] = ""
        return df

    # Load taxonomy
    taxonomy_df = pd.read_csv(
        taxonomy,
        sep="\t",
        dtype=str,
        comment="#",
        low_memory=False,
    )

    # Clean column names
    taxonomy_df.columns = (
        taxonomy_df.columns
        .str.strip()
        .str.lower()
    )

    # Fix known header typo
    taxonomy_df = taxonomy_df.rename(
        columns={"phylumphylum_taxon_id": "phylum_taxon_id"}
    )

    # Validate required columns
    if "accession" not in taxonomy_df.columns or "family" not in taxonomy_df.columns:
        raise ValueError(
            f"Missing required columns. Found: {taxonomy_df.columns.tolist()}"
        )

    # Build lookup table
    lookup_df = taxonomy_df[["accession", "family"]].copy()

    lookup_df["accession"] = (
        lookup_df["accession"]
        .fillna("")
        .astype(str)
        .str.strip()
        .str.upper()
        .str.replace(r"\.\d+$", "", regex=True)  # remove version suffix
    )

    lookup_df["family"] = (
        lookup_df["family"]
        .fillna("")
        .astype(str)
        .str.strip()
    )

    # Drop duplicates after normalisation
    lookup_df = lookup_df.drop_duplicates(subset=["accession"], keep="first")

    accession_map = dict(zip(lookup_df["accession"], lookup_df["family"]))

    # Apply to input df
    enriched_df = df.copy()

    enriched_df["accession"] = (
        enriched_df["accession"]
        .fillna("")
        .astype(str)
        .str.strip()
        .str.upper()
        .str.replace(r"\.\d+$", "", regex=True)  # ensure same format
    )

    enriched_df["family"] = (
        enriched_df["accession"]
        .map(accession_map)
        .fillna("")
    )

    return enriched_df


def filter_support(df, min_reads):
    work = df.copy()
    work["reads"] = pd.to_numeric(work.get("reads"), errors="coerce").fillna(0).astype(int)
    target_levels = ["species_unclassified", "genus_unclassified"]
    resolution = work.get("resolution_level", pd.Series([""] * len(work)))
    return work[
        (resolution.isin(target_levels)) & (work["reads"] > min_reads)
    ].copy()



def load_diamond_results(path):
    """Load diamond results based on mode."""
    #if not file_exists(path):
    #    raise FileNotFoundError(f"Diamond results file not found: {path}")
     # Define expected header (based on your Diamond output fields)
    columns = [
        "qseqid", "sseqid", "pident", "alignment_length", "qlen", "slen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "stitle"
    ]
        # Create a temporary file with header + original content
    df = pd.read_csv(path, sep="\t", header=None, dtype=str)
    if df.shape[1] != len(columns):
        raise ValueError(
            f"Expected {len(columns)} columns (BLAST output), but found {df.shape[1]} in {path}"
        )
    df.columns = columns

    dtype = {
        "qseqid": 'str', "sseqid": 'str', "pident": 'float64', "alignment_length": 'int64', 
        "qlen": 'int64', "slen": 'int64', "qstart": 'int64', "qend": 'int64', "sstart": 'int64', 
        "send": 'int64', "evalue": 'float64', "bitscore": 'int64', "stitle": 'str'
    }

    # Load DataFrame
    for col, dtype_ in dtype.items():
        if col in df.columns:
            if dtype_ in ("int64", "float64"):
                df[col] = pd.to_numeric(df[col], errors="coerce").astype(dtype_)
            else:
                df[col] = df[col].astype(str)
    df = df[df["stitle"].str.contains("virus", case=False, na=False)].copy()
    df = df[
        ~df["stitle"].str.contains("retrovirus-related|mitovirus|hypothetical|putative|retrovirus|Caulimovirus sp|mito-like|ambiguivirus|narna-like|picorna-like|Plasmopara|Mimivirus|tombus-like|Solanum tuberosum", case=False, na=False)
    ].copy()
    #df["stitle_short"] = df["stitle"].astype(str).str.rsplit("|", n=1).str[-1].str.strip()
    df["taxon"] = (
        df["stitle"]
        .astype(str)
        .str.extract(r"\[([^\[\]]+)\]\s*$")[0]   # last [...] at end of string
        .fillna("")
        .str.strip()
    )
    df["desc"] = (
        df["stitle"]
        .astype(str)
        .str.split("|")
        .str[-1]                              # last pipe element only
        .str.replace(r"\s*\[[^\[\]]+\]\s*$", "", regex=True)  # remove trailing [..]
        .str.strip()
    )
    df["accession"] = (
        df["stitle"]
        .astype(str)
        .str.split("|")
        .str[-2]
        .fillna("")
        .str.strip()
    )
    
    df["qcovs"] = (df["alignment_length"] * 3  / df["qlen"]) * 100
    #df = (df["qcovs"] > 75).astype(int)
    df = df[df["qcovs"] >= 50].copy()
    df = df[df["pident"] >= 50].copy()
    df = df.sort_values(by="pident", ascending=False)
    df2 = df[(df["pident"] < 90) & (df["qlen"] >= 500)].copy()

    final_columns_filt = [
        "qseqid",
        "taxon",
        "accession",
        "desc",
        "pident",
        "qlen",
        "alignment_length",
        "evalue",
        "bitscore",
        "qcovs",
    ]
    
    return df[final_columns_filt], df2[final_columns_filt]

def build_rows(Method, filtered_df):
    if filtered_df.empty:
        return [
            {
                "Method": Method,
                "Evidence": "",
                "Details": "0 reads",
                "Taxonomy_classification": "",
            }
        ]

    grouped = (
        filtered_df.assign(
            taxon_name=filtered_df.get("taxon_name", "").fillna("").astype(str).str.strip(),
            resolution_level=filtered_df.get("resolution_level", "").fillna("").astype(str).str.strip(),
            reads=filtered_df["reads"].astype(int),
        )
        .groupby(["resolution_level", "taxon_name"], dropna=False, as_index=False)["reads"]
        .sum()
    )

    level_labels = {
        "genus_unclassified": "genus unresolved",
        "species_unclassified": "species unresolved",
    }

    rows = []
    for _, row in grouped.iterrows():
        taxon = row["taxon_name"]
        reads = int(row["reads"])
        level_label = level_labels.get(row["resolution_level"], row["resolution_level"])

        rows.append(
            {
                "Method": Method,
                "Evidence": "High read counts assigned to viral taxa at higher taxonomic level",
                "Details": f"{reads} reads",
                "Taxonomy_classification": f"{taxon} ({level_label})",
            }
        )


    return rows

def build_megablast_rows(megablast_df):
    work = megablast_df.copy()
    work["pident"] = pd.to_numeric(work.get("pident"), errors="coerce").fillna(0)
    work["qlen"] = pd.to_numeric(work.get("qlen"), errors="coerce").fillna(0)
    work["mapping_read_count"] = pd.to_numeric(
        work.get("mapping_read_count"), errors="coerce"
    ).fillna(0).astype(int)
    work["taxon_name"] = work.get("taxon_name", "").fillna("").astype(str).str.strip()
    work["qseqid"] = work.get("qseqid", "").fillna("").astype(str).str.strip()
    work["qlen"] = pd.to_numeric(work.get("qlen"), errors="coerce").fillna(0).astype(int)
    work["pc_cov_30X"] = pd.to_numeric(work.get("pc_cov_30X"), errors="coerce").fillna(0).astype(float)
    work["normalised_conf_score"] = pd.to_numeric(work.get("normalised_conf_score"), errors="coerce").fillna(0).astype(float)
    #extract lineage from family downwards
    work["lineage"] = (
        work["full_lineage"]
        .fillna("")
        .astype(str)
        .apply(lambda x: ";".join(x.split(";")[6:]) if ";" in x else "")
    )
   
    filtered = work[(work["pident"] <= 90) & (work["qlen"] >= 1000) & (work["pc_cov_30X"] >= 70)].copy()

    if filtered.empty:
        return [
            {
                "Method": "Megablast",
                "Evidence": "",
                "Details": "",
                "Taxonomy_classification": "",
            }
        ]

    grouped = (
        filtered.assign(pident=filtered["pident"].round(2))
        .groupby(["lineage", "taxon_name", "pident", "qseqid", "qlen", "pc_cov_30X", "normalised_conf_score"], as_index=False)["mapping_read_count"]
        .sum()
    )

    rows = []
    for _, row in grouped.iterrows():
        rows.append(
            {
                "Method": "Megablast",
                "Evidence": f"Support for low-identity viral match",
                "Details": f"{row['qseqid']}; {row['pident']}% id; {row['qlen']} nt; {int(row['mapping_read_count'])} mapping reads; {row['pc_cov_30X']}% 30X; {row['normalised_conf_score']} norm_conf_score",
                "Taxonomy_classification": f"{row['lineage']}",
            }
        )

    return rows
#"Details": f"contig={row['qseqid']}; pident={row['pident']}, qlen={row['qlen']}, reads={int(row['mapping_read_count'])}, pc_cov_30X={row['pc_cov_30X']}, normalised_conf_score={row['normalised_conf_score']}",

def build_diamond_rows(diamond_df):
    if diamond_df.empty:
        return [
            {
                "Method": "Diamond",
                "Evidence": "",
                "Details": "",
                "Taxonomy_classification": "",
            }
        ]

    rows = []
    for _, row in diamond_df.iterrows():
        #family = str(row.get("family", "")).strip()
        #taxon = str(row.get("taxon", "")).strip()
        #taxonomy_label = f"{family}; {taxon}" if family else taxon
        rows.append(
            {
                "Method": "Diamond",
                "Evidence": f"Low-identity viral match with contig length >= 500 nt",
                "Details": f"{row['qseqid']}; {row['pident']}% id; {row['qlen']} nt; {row['qcovs']:.1f}% qcovs",
                "Taxonomy_classification": f"{row['family']}; {row['taxon']}",
            }
        )
    return rows


def build_novel_rows(novel_df):
    work = novel_df.copy()
    work["virus_score"] = pd.to_numeric(work.get("virus_score"), errors="coerce").fillna(0)
    work["ORFs"] = pd.to_numeric(work.get("ORFs"), errors="coerce").fillna(0)
    work["RdRp"] = pd.to_numeric(work.get("RdRp"), errors="coerce").fillna(0)
    work["PFAM_total"] = pd.to_numeric(work.get("PFAM_total"), errors="coerce").fillna(0)
    work["taxonomy"] = work.get("taxonomy", "").fillna("").astype(str).str.strip()

    def taxonomy_last(value):
        parts = [x.strip() for x in str(value).split(";") if x and x.strip()]
        return parts[-1] if parts else ""

    filtered = work[
        (work["virus_score"] > 0.99) & (work["ORFs"] >= 1) & (work["RdRp"] >= 1)
    ].copy()

    if filtered.empty:
        return [
            {
                "Method": "Genomad",
                "Evidence": "",
                "Details": "",
                "Taxonomy_classification": "",
            }
        ]

    rows = []
    for _, row in filtered.iterrows():
        rows.append(
            {
                "Method": "Genomad",
                "Evidence": "Functional evidence for a putative novel virus in a contig returning no megablast hit",
                "Details": (
                    f"{row['seq_name']}; "
                    f"virus_score={row['virus_score']:.3f}; "
                    f"ORFs={int(row['ORFs'])}; RdRp={int(row['RdRp'])}; "
                    f"PFAM_total={int(row['PFAM_total'])}"
                ),
                "Taxonomy_classification": taxonomy_last(row["taxonomy"]),
            }
        )
    return rows



def main():
    ################################################################################
    parser = argparse.ArgumentParser(description="Load blast and coverage stats summary")
    parser.add_argument("--sample_name", type=str, required=True, help='provide sample name')
    parser.add_argument("--blast", type=str, required=True, help='provide blast top hits')
    parser.add_argument("--kraken", type=str, required=True, help='provide kraken file')
    parser.add_argument("--kaiju", type=str, required=True, help='provide kaiju file')
    parser.add_argument("--hmmscan", type=str, required=True, help='provide hmmscan file')
    parser.add_argument("--map2ref", type=str, required=True, help='provide coverage stats to reference')
    parser.add_argument("--fasta", type=str, required=True, help='provide fasta file')
    parser.add_argument("--genomad", type=str, required=True, help='provide genomad results')
    parser.add_argument("--blast_novel", type=str, required=False, help='provide output file name')
    parser.add_argument("--diamond", required=True, help="Diamond matches file (*diamond_matches.txt)")
    parser.add_argument("--min-reads", type=int, default=2000, help="Strict minimum reads threshold (reads > min_reads)")
    parser.add_argument("--taxonomy", required=True, help="Path to the taxonomy file")

    args = parser.parse_args()
    sample_name = args.sample_name
    blast = args.blast
    kraken = args.kraken
    kaiju = args.kaiju
    hmmscan = args.hmmscan
    map2ref = args.map2ref
    fasta = args.fasta
    genomad = args.genomad
    blast_novel = args.blast_novel
    diamond = args.diamond
    min_reads = args.min_reads
    taxonomy = args.taxonomy
    

    df = parse_hmmscan_per_target(hmmscan)
    df_filtered = filter_hmmscan(df)
    print("Number of matching hits:", len(df_filtered))

    hmm_df = summarize_domains(df_filtered)
    hmm_df.to_csv(sample_name + "_hmm_domain_summary_counts.tsv", sep="\t", index=False)


    blast_df = pd.read_csv(blast, sep="\t", dtype=str)
    kraken_df = pd.read_csv(kraken, sep="\t", dtype=str)
    kaiju_df = pd.read_csv(kaiju, sep="\t", dtype=str)
    map2ref_df = pd.read_csv(map2ref, sep="\t", dtype=str)
    print(blast_df.head())

    #filter blast results to only those that passed all filters inc best contig per species filter
    blast_df["pident"] = pd.to_numeric(blast_df["pident"], errors="coerce")
    blast_df["qlen"] = pd.to_numeric(blast_df["qlen"], errors="coerce")
    blast_df["mapping_read_count"] = pd.to_numeric(blast_df["mapping_read_count"], errors="coerce")
    blast_df["pc_mapping_reads"] = pd.to_numeric(blast_df["pc_mapping_reads"], errors="coerce")
    blast_df = add_group_max(blast_df, "species", "pident", "max_pident_spp")
    blast_df = add_group_max(blast_df, "species", "qlen", "max_qlen_spp")
    blast_df = add_group_max(blast_df, "species", "mapping_read_count", "max_contig_mapping_read_count_spp")
    blast_df = add_group_max(blast_df, "species", "pc_mapping_reads", "max_pc_contig_mapping_reads_spp")


    #only reports PFAMs for contigs that passed all filters and were the best contig per species
    #chnage to stre whether there was a PFAM hit for any contig per species that passed filters, not just the best contig per species
    blast_df2 = blast_df.merge(
        hmm_df[["query_name", "PFAM_total", "RdRp"]],
        left_on="qseqid",
        right_on="query_name",
        how="left"   # left join keeps all rows in summary_df
    )
    blast_df2 = add_group_max(blast_df2, "species", "PFAM_total", "max_PFAM_total_spp")
    # ensure RdRp is boolean
    blast_df2["RdRp"] = blast_df2["RdRp"].fillna(False).infer_objects(copy=False)

    # species-level flag: any contig of this species has RdRp
    blast_df2["species_has_RdRp"] = (
        blast_df2.groupby("species")["RdRp"]
        .transform("any")
    )
    

    filtered_blast_df = blast_df2[
        (blast_df2["term_filter"].astype(str) == "True") &
        (blast_df2["cov_filter"].astype(str) == "True") &
        (blast_df2["best_contig_per_sp_filter"].astype(str) == "True")
    ]

    # Standardize species naming across all tables for the same viroid
    filtered_blast_df = standardize_viroid_names(filtered_blast_df, "species", df_name="filtered_blast_df")
    kraken_df = standardize_viroid_names(kraken_df, "taxon_name", df_name="kraken_df")
    kaiju_df = standardize_viroid_names(kaiju_df, "taxon_name", df_name="kaiju_df")
    map2ref_df = standardize_viroid_names(map2ref_df, "taxon_name", df_name="map2ref_df")

    #print(kraken_df.head())
    #print(kaiju_df.head())
    summary_df = pd.DataFrame()
    
    
    # Mask Kraken/Kaiju species to only viral entries and filter for those that passed the term filter
    kraken_viral = kraken_df.loc[
        (kraken_df["broad_category"] == "viral") & (kraken_df["term_filter"].str.lower() == "true"),
        "taxon_name"
    ]
    kaiju_viral  = kaiju_df.loc[
        (kaiju_df["broad_category"] == "viral") & (kaiju_df["term_filter"].str.lower() == "true"),
        "taxon_name"
    ]

    # Combine vertically and get unique
    summary_df["taxon"] = pd.Series(
        pd.unique(
            pd.concat([
                filtered_blast_df["species"],
                kraken_viral,
                kaiju_viral
            ], ignore_index=True)
        )
    )
    # Add boolean column: True if species appeared in blast
    blast_species_set = set(filtered_blast_df["species"])
    kraken_species_set = set(kraken_df["taxon_name"])
    kaiju_species_set = set(kaiju_df["taxon_name"])
    summary_df["taxon_in_blast"] = summary_df["taxon"].isin(blast_species_set)
    summary_df["taxon_in_kraken"] = summary_df["taxon"].isin(kraken_species_set)
    summary_df["taxon_in_kaiju"] = summary_df["taxon"].isin(kaiju_species_set)
    #print(summary_df)

    merged_df = summary_df.merge(
        filtered_blast_df[["species", "qseqid", "qlen", "sacc", "pident", "bitscore", "evalue", "contig_seq", "ncontigs_per_spp", 
                           "max_pident_spp","max_qlen_spp", "max_contig_mapping_read_count_spp",  "max_pc_contig_mapping_reads_spp", "total_score_spp", "mapping_read_count", "pc_mapping_reads", "mean_depth", "pc_cov_30X", 
                            "mean_mapping_quality", "read_count_flag", "mean_depth_flag", "30x_cov_flag", "mean_mq_flag", "species_has_RdRp", "max_PFAM_total_spp",
                            "total_conf_score","normalised_conf_score"]],
        left_on="taxon",
        right_on="species",
        how="left"   # left join keeps all rows in summary_df
    )
    merged_df.rename(columns={"qlen": "contig_length",
                               "contig_seq": "contig_seq"}, inplace=True)


    merged_df2 = merged_df.merge(
        kaiju_df[["taxon_name", "reads", "pc_reads"]],
        left_on="taxon",
        right_on="taxon_name",
        how="left"   # left join keeps all rows in summary_df
    )
    merged_df2[["reads", "pc_reads"]] = (
        merged_df2[["reads", "pc_reads"]].fillna(0)
    )
    # Rename the columns
    merged_df2.rename(columns={"pc_reads": "kaiju_pc_reads",
                               "reads": "kaiju_reads"}, inplace=True)


    merged_df3 = merged_df2.merge(
        kraken_df[["taxon_name", "reads", "pc_reads"]],
        left_on="taxon",
        right_on="taxon_name",
        how="left"   # left join keeps all rows in summary_df
    )
    merged_df3[["reads", "pc_reads"]] = (
        merged_df3[["reads", "pc_reads"]].fillna(0)
    )
    # Rename the columns
    merged_df3.rename(columns={"pc_reads": "kraken_pc_reads",
                               "reads": "kraken_reads"}, inplace=True)

    #ow many times does this taxon_name and sacc pairing appear
    #map2ref_df["ref_count"] = map2ref_df.groupby(["taxon_name", "sacc"])["sacc"].transform("count")
    map2ref_df["mapping_read_count"] = pd.to_numeric(
        map2ref_df["mapping_read_count"], errors="coerce"
    )
    map2ref_df = add_group_max(map2ref_df, "taxon_name", "mapping_read_count", "max_ref_mapping_read_count_spp")
    map2ref_df = add_group_max(map2ref_df, "taxon_name", "pc_mapping_reads", "max_pc_ref_mapping_reads_spp")
    map2ref_df = add_group_max(map2ref_df, "taxon_name", "pc_cov_30X", "max_pc_cov_30X_spp")
    map2ref_df = add_group_max(map2ref_df, "taxon_name", "normalised_conf_score", "max_normalised_conf_score_spp")
    map2ref_df = add_group_max(map2ref_df, "taxon_name", "reference_length", "max_reference_length_spp")
    map2ref_df_unique = (
        map2ref_df[
            ["taxon_name",
            "max_ref_mapping_read_count_spp",
            "max_pc_ref_mapping_reads_spp" ,
            "max_reference_length_spp",
            "max_normalised_conf_score_spp",
            "max_pc_cov_30X_spp"]
            ]
            .drop_duplicates(subset=["taxon_name"])
        )
    
    merged_df4 = merged_df3.merge(
        map2ref_df_unique[["taxon_name", "max_ref_mapping_read_count_spp", "max_pc_ref_mapping_reads_spp", "max_reference_length_spp", "max_normalised_conf_score_spp", "max_pc_cov_30X_spp"]],
        left_on="taxon",
        right_on="taxon_name",
        how="left"   # left join keeps all rows in summary_df
    )

    final_columns_filt = ["taxon","kraken_reads", "kraken_pc_reads", 
                          "max_contig_mapping_read_count_spp",  "max_pc_contig_mapping_reads_spp",
                          "max_ref_mapping_read_count_spp", "max_pc_ref_mapping_reads_spp",
                          "kaiju_reads", "kaiju_pc_reads", 
                          "mapping_read_count","pc_mapping_reads",
                          "ncontigs_per_spp", "qseqid", "contig_seq", "max_pident_spp", "max_qlen_spp",
                          "mean_depth", "pc_cov_30X", "normalised_conf_score", "max_PFAM_total_spp", "species_has_RdRp", 
                          "max_reference_length_spp", "max_normalised_conf_score_spp", "max_pc_cov_30X_spp"]

    merged_df4 = merged_df4[final_columns_filt]
    text_cols = ["qseqid", "contig_seq","species_has_RdRp" ]
    #text_cols = ["qseqid", "contig_seq", "sacc", "read_count_flag", "mean_depth_flag", "30x_cov_flag", "mean_mq_flag", "RdRp" ]
    num_cols = ["ncontigs_per_spp","mapping_read_count", "pc_mapping_reads", "mean_depth", "pc_cov_30X", "normalised_conf_score" , "max_PFAM_total_spp",  "max_pident_spp", "max_qlen_spp"]

    # fill values
    merged_df4 = merged_df4[merged_df4["contig_seq"].notna()]
    merged_df4[text_cols] = merged_df4[text_cols].fillna("NA")
    merged_df4[num_cols] = merged_df4[num_cols].fillna(0)
    
    merged_df4[num_cols] = merged_df4[num_cols].apply(pd.to_numeric, errors="coerce")

    merged_df4.sort_values(
        by=["mapping_read_count", "normalised_conf_score"],
        ascending=[False, False],
        inplace=True
    )
    
    output_file = f"{sample_name}_summary_viral_results.tsv"
    merged_df4.to_csv(output_file, index=False, sep="\t")

    #derive novel virus candidates summary table, by merging the fasta, genomad, hmmscan and blast results and applying filters to identify novel candidates
    if file_exists(fasta) and not is_file_empty(fasta):
        fasta_df = fasta_to_dataframe(fasta)
    else:
        fasta_df = pd.DataFrame(columns=["seq_name", "contig_seq"])

    if file_exists(genomad) and not is_file_empty(genomad):
        
        genomad_results = pd.read_csv(genomad, sep="\t", header=0)
        fasta_genomad_df = pd.merge(fasta_df, genomad_results, on = ['seq_name'], how = 'outer')


    #if file_exists(hmmscan) and not is_file_empty(hmmscan):
        
    #    hmmscan_results = pd.read_csv(hmmscan, sep="\t", header=0)
    fasta_genomad_hmmscan_df = pd.merge(fasta_genomad_df, hmm_df, 
            left_on = ['seq_name'], 
            right_on = ['query_name'],
            how = 'left')
        
    if file_exists(blast_novel) and not is_file_empty(blast_novel):
        
        blast_results = pd.read_csv(blast_novel, sep="\t", header=0)
        fasta_genomad_hmmscan_blast_df = pd.merge(fasta_genomad_hmmscan_df, blast_results, 
                              left_on = ['seq_name'], 
                              right_on = ['qseqid'],
                              how = 'outer')

    num_cols = ["length", "virus_score", "PFAM_total"]
    fasta_genomad_hmmscan_blast_df[num_cols] = fasta_genomad_hmmscan_blast_df[num_cols].apply(pd.to_numeric, errors="coerce")
    # Force integer columns
    int_cols = ["n_genes", "length", "n_hallmarks", "genetic_code"]
    for col in int_cols:
        fasta_genomad_hmmscan_blast_df[col] = pd.to_numeric(fasta_genomad_hmmscan_blast_df[col], errors="coerce").astype("Int64")

    rows_to_check = ["PFAM_total", "virus_score"]

    df_filt = fasta_genomad_hmmscan_blast_df[
        ~(fasta_genomad_hmmscan_blast_df[rows_to_check].fillna(0).eq(0).all(axis=1))
    ]
    #consider filtering the unclassified hits, as these are likely to be false positives?
    df_filt = df_filt[
        ~df_filt["taxonomy"].str.lower().str.contains("caudoviricetes|iridoviridae|retroviridae", na=False)
    ]
    df_filt = df_filt[df_filt["length"] >= 500]
    df_filt = df_filt[df_filt["virus_score"] >= 0.99]
    #Keep only contigs with no blast hits.
    df_filt = df_filt[
        df_filt["sacc"].isna() | (df_filt["sacc"].str.strip() == "")
    ]


    df_filt["ORFs"] = df_filt["ORFs"].fillna("NA")

    # Convert RdRp to numeric (True=1, False=0)
    df_filt["RdRp"] = (
        df_filt["RdRp"]
        .fillna(False)
        .astype(bool)
        .astype(int)
    )

    # Drop low-confidence unclassified contigs lacking ORFs and RdRp support.
    orfs_numeric = pd.to_numeric(df_filt["ORFs"], errors="coerce").fillna(0)
    df_filt = df_filt[
        ~(
            df_filt["taxonomy"].str.strip().str.lower().eq("unclassified")
            & orfs_numeric.eq(0)
            & df_filt["RdRp"].eq(0)
        )
    ]

    # FORCE correct dtypes before sorting
    df_filt["virus_score"] = pd.to_numeric(
        df_filt["virus_score"], errors="coerce"
    ).fillna(0)

    df_filt["PFAM_total"] = pd.to_numeric(
        df_filt["PFAM_total"], errors="coerce"
    ).fillna(0)


    final_columns_filt = ["seq_name", "contig_seq", "length", "n_genes", "genetic_code", "virus_score", "n_hallmarks", "marker_enrichment", "taxonomy", "ORFs", "RdRp", "PFAM_total"]

    df_filt = df_filt[final_columns_filt]
    df_filt.sort_values(
        by=["virus_score", "RdRp", "PFAM_total"],
        ascending=[False, False, False],
        inplace=True
    )

    output_file = f"{sample_name}_novel_virus_candidates.tsv"
    df_filt.to_csv(output_file, index=False, sep="\t")

    #derive summary table for evidence of novel viruses, by merging the kaiju, kraken, megablast, hmmscan and diamond results and applying filters to identify novel candidates. This is more focused on summarising evidence for novel viruses, whereas the previous table is focused on summarising all viral hits including known viruses.
    kaiju_df = pd.read_csv(kaiju, sep="\t", dtype=str)
    kraken_df = pd.read_csv(kraken, sep="\t", dtype=str)
    megablast_df = pd.read_csv(blast, sep="\t", dtype=str)
    #novel_df = pd.read_csv(args., sep="\t", dtype=str)
    #diamond_df = pd.read_csv(args.diamond, sep="\t", dtype=str)
    diamond_results_path = diamond
    diamond_results, diamond_novel_candidate_results = load_diamond_results(diamond_results_path)
    diamond_results.to_csv(f"{args.sample_name}_filtered_diamond_results.txt", sep="\t", index=False)
    #diamond_novel_candidate_results.to_csv(f"{args.sample_name}_novel_candidate_diamond_results.txt", sep="\t", index=False)
    kaiju_filtered = filter_support(kaiju_df, min_reads)
    kraken_filtered = filter_support(kraken_df, min_reads)

    support_rows = []
    support_rows.extend(build_rows("Kaiju", kaiju_filtered))
    support_rows.extend(build_rows("Kraken", kraken_filtered))
    support_rows.extend(build_megablast_rows(megablast_df))
    enriched_df = enrich_with_taxonomy(diamond_novel_candidate_results, taxonomy)
    #print(enriched_df)
    support_rows.extend(build_diamond_rows(enriched_df))
    support_rows.extend(build_novel_rows(df_filt))
    support_df = pd.DataFrame(support_rows)
    support_df.to_csv(f"{args.sample_name}_evidence_summary_novel.txt", sep="\t", index=False)

if __name__ == "__main__":
    main()