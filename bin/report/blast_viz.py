#!/usr/bin/env python

import base64
import io
import re
import os
import html
import pandas as pd


# --- Utility: BLAST Parsing to get all reference & per-contig alignment boxes (one per contig/ref) ---


def svg_to_b64(svg_string):
    return base64.b64encode(svg_string.encode()).decode()

def parse_blast_for_ruler_continuous(blast_filename):
    """
    Parses BLAST output and returns:
        - ref_lengths: {ref_id: length}
        - merged_contig_alignments: [{'contig': ..., 'ref': ..., 'ref_start': ..., 'ref_end': ...}]
    Merges all blocks for each (contig, ref) into a single interval (so one box per contig/ref).
    """
    ref_lengths = {}
    contig_blocks = {}  # key = (contig, ref), value = [all starts, all ends]
    current_ref = None
    current_contig = None

    with open(blast_filename) as f:
        lines = f.readlines()
    n = len(lines)
    i = 0
    while i < n:
        line = lines[i]
        m_q = re.match(r"^\s*Query=\s*(\S+)", line)
        if m_q:
            current_contig = m_q.group(1)
        m_ref = re.match(r"^>(\S+)", line)
        if m_ref:
            current_ref = m_ref.group(1)
            k = i+1
            while k < n and not lines[k].startswith(">") and lines[k].strip() != "":
                m_len = re.match(r"Length\s*=\s*(\d+)", lines[k])
                if m_len:
                    ref_clean = current_ref.split('.')[0]
                    ref_lengths[ref_clean] = int(m_len.group(1))
                    break
                k += 1
        if line.startswith("Sbjct"):
            m_s = re.match(r"Sbjct\s+(\d+)\s+[A-Za-z\-]+\s+(\d+)", line)
            if m_s and current_ref and current_contig:
                start, end = int(m_s.group(1)), int(m_s.group(2))
                key = (current_contig, current_ref)
                if key not in contig_blocks:
                    contig_blocks[key] = {"starts": [], "ends": [], "contig": current_contig, "ref": current_ref}
                contig_blocks[key]["starts"].append(start)
                contig_blocks[key]["ends"].append(end)
        i += 1
    merged_contig_alignments = []
    for v in contig_blocks.values():
        all_pos = v["starts"] + v["ends"]
        merged_contig_alignments.append({
            "contig": v["contig"],
            "ref": v["ref"],
            "ref_start": min(all_pos),
            "ref_end":   max(all_pos)
        })
    return ref_lengths, merged_contig_alignments

# --- Utility: VirusDetect-style table with alignments ---
def extract_top_hit_table_with_alignment(blast_file):
    results = []
    with open(blast_file, 'r') as f:
        lines = f.readlines()
    n = len(lines)
    i = 0
    while i < n:
        m_query = re.match(r"\s*Query=\s*(\S+)", lines[i])
        if m_query:
            query = m_query.group(1)
            hit = None; score = evalue = identity = hsp_length = None
            query_start = query_end = hit_start = hit_end = ""
            strand = ""; identity2 = ""
            aligned_query = []; aligned_hit = []; aligned_string = []; alignment_lines = []
            # 1. Find top hit from "Sequences producing significant alignments:"
            i2 = i+1; in_significant = False
            while i2 < n:
                if "Sequences producing significant alignments:" in lines[i2]:
                    in_significant = True; i2 += 1; continue
                if in_significant:
                    if lines[i2].strip() and not lines[i2].startswith(">"):
                        hit = lines[i2].split()[0]; break
                i2 += 1
            if not hit: i += 1; continue
            # 2. Find >hit section
            i3 = i2
            while i3 < n:
                if lines[i3].startswith(">" + hit): break
                i3 += 1
            # 3. Find first Score (top HSP)
            i4 = i3
            while i4 < n:
                m_score = re.match(r".*Score =\s*([\d\.]+) bits.*Expect.* = (\S+)", lines[i4])
                if m_score:
                    score = m_score.group(1); evalue = m_score.group(2)
                    j = i4 + 1
                    while j < n:
                        m_id = re.match(r"\s*Identities = (\d+)/(\d+) \(([\d\.]+)%\).*", lines[j])
                        if m_id:
                            identity = m_id.group(3); hsp_length = m_id.group(2)
                            identity2 = f"{m_id.group(1)}/{m_id.group(2)}({m_id.group(3)}%)"
                        m_strand = re.match(r"\s*Strand=([^\n]*)", lines[j])
                        if m_strand:
                            strand = m_strand.group(1).strip()
                        if re.match(r"^Query\s", lines[j]): break
                        j += 1
                    q_min = q_max = h_min = h_max = None; k = j
                    while k < n:
                        l = lines[k]
                        m_q = re.match(r"^Query\s+(\d+)\s+([A-Za-z\-]+)\s+(\d+)", l)
                        if m_q:
                            qs, qseq, qe = int(m_q.group(1)), m_q.group(2), int(m_q.group(3))
                            if q_min is None or qs < q_min: q_min = qs
                            if q_max is None or qe > q_max: q_max = qe
                            aligned_query.append(qseq)
                            # collect match line and Sbjct line if present
                            if k+1 < n:
                                aligned_string.append(lines[k+1].rstrip())
                            if k+2 < n:
                                m_s = re.match(r"^Sbjct\s+(\d+)\s+([A-Za-z\-]+)\s+(\d+)", lines[k+2])
                                if m_s:
                                    hs, hseq, he = int(m_s.group(1)), m_s.group(2), int(m_s.group(3))
                                    if h_min is None or hs < h_min: h_min = hs
                                    if h_max is None or he > h_max: h_max = he
                                    aligned_hit.append(hseq)
                                else:
                                    aligned_hit.append("")
                            else:
                                aligned_hit.append("")
                            # Collect lines for the HTML block
                            alignment_lines.append(l.rstrip())
                            if k+1 < n: alignment_lines.append(lines[k+1].rstrip())
                            if k+2 < n: alignment_lines.append(lines[k+2].rstrip())
                            alignment_lines.append("")
                            k += 3
                            continue
                        # Stop at new hit, query, or HSP
                        if (l.startswith(">") or l.startswith("Query=") or 
                            "Sequences producing significant alignments:" in l or 
                            re.match(r".*Score =\s*[\d\.]+ bits", l)):
                            break
                        k += 1

                    results.append([
                        query, hit, score if score else "", evalue if evalue else "", identity if identity else "",
                        hsp_length if hsp_length else "", str(q_min) if q_min is not None else "",
                        str(q_max) if q_max is not None else "", str(h_min) if h_min is not None else "",
                        str(h_max) if h_max is not None else "", strand, identity2,
                        ",".join(aligned_query), ",".join(aligned_hit), ",".join(aligned_string),
                        "\n".join(alignment_lines).rstrip("\n")
                    ])
                    break
                i4 += 1
        i += 1
    return results


# --- HTML Table Renderer (VirusDetect style) ---
def output_alignments_html(align_info, svg_filename=None):
    """
    align_info: rows from extract_top_hit_table_with_alignment
      [
        query, hit, score, evalue, identity, hsp_length,
        query_start, query_end, hit_start, hit_end, strand,
        identity2, aligned_query, aligned_hit, aligned_string, alignment_block
      ]
    """
    out_html = ""
    if svg_filename:
        out_html += f'''
<div style="text-align:center; margin-bottom:18px;">
  <img src="{svg_filename}" alt="Reference Plot" style="max-width:98%;">
</div>
'''
    out_html += '''
<table border="1" cellspacing="0" cellpadding="6" style="border-collapse:collapse; font-family:sans-serif; font-size:13px; width:100%;">
  <tr bgcolor="#e2e8ec">
    <th>Order</th>
    <th>Query ID</th>
    <th>Query Start</th>
    <th>Query End</th>
    <th>Subject Start</th>
    <th>Subject End</th>
    <th>Identity</th>
    <th>E value</th>
    <th>Strand</th>
  </tr>
'''

    for idx, row in enumerate(align_info, 1):
        query, hit, score, evalue, identity, hsp_length, qs, qe, hs, he, strand, identity2, aligned_query, aligned_hit, aligned_string, alignment_block = row
        # Table row
        out_html += f'''
<tr>
  <td>{idx}</td>
  <td style="color:#178cd2;"><a id="{query}">{query}</a></td>
  <td>{qs}</td>
  <td>{qe}</td>
  <td>{hs}</td>
  <td>{he}</td>
  <td>{identity2}</td>
  <td>{evalue}</td>
  <td>{strand}</td>
</tr>
<tr>
  <td colspan="9" style="background:#fcfcfc; border-top:0;">
    <div>
      <b>Alignment:</b>
      <pre style="background:#f0fdff; color:#333; padding:8px; border-radius:4px; font-size:13px;">
{highlight_alignment_block(alignment_block)}
      </pre>
    </div>
  </td>
</tr>
'''
    out_html += '</table>'
    return out_html

def highlight_alignment_block(alignment_block):
    """
    HTML-color the full multi-line alignment block.
    Query and Sbjct lines get background highlight, as in VirusDetect.
    """
    out = []
    for line in alignment_block.rstrip('\n').split('\n'):
        line = html.escape(line)
        if line.startswith("Query"):
            out.append(f'<span style="background-color:#c4fcfc;">{line}</span>')
        elif line.startswith("Sbjct"):
            out.append(f'<span style="background-color:#c4fcfc;">{line}</span>')
        elif set(line.strip()) <= set("| "):  # match line (all bars and spaces)
            out.append(line)
        elif line.strip() == "":
            out.append("")
        else:
            out.append(line)
    return "\n".join(out)


def generate_reference_svg(ref_length, contigs):
    width = 900   # or 1000 for safety
    right_padding = 100
    scale = (width - right_padding) / ref_length
    height = 50 + len(contigs) * 20

    scale = width / ref_length

    svg = [f'<svg width="{width}" height="{height}" xmlns="http://www.w3.org/2000/svg">']

    # Reference line
    ref_px = ref_length * scale

    svg.append(
        f'<line x1="0" y1="20" x2="{ref_px}" y2="20" stroke="#1267b2" stroke-width="14"/>'
    )

    # Ticks
    for pos in range(0, ref_length + 1, 1000):
        x = pos * scale
        svg.append(f'<line x1="{x}" y1="15" x2="{x}" y2="25" stroke="black"/>')
        svg.append(f'<text x="{x}" y="40" font-size="10" text-anchor="middle">{pos//1000}k</text>')

    # Contigs
    for i, c in enumerate(contigs):
        y = 60 + i * 20
        x1 = c['ref_start'] * scale
        x2 = c['ref_end'] * scale

        svg.append(
            f'<rect x="{x1}" y="{y}" width="{x2-x1}" height="16" fill="#e33"/>'
        )
        x_center = (x1 + x2) / 2

        # If label would overflow right edge
        if x_center > width - 60:
            x_text = x2            # anchor to end of contig
            anchor = "end"
        else:
            x_text = x_center
            anchor = "middle"

        svg.append(
            f'<text x="{x_text}" y="{y-2}" font-size="9" text-anchor="{anchor}">{c["contig"]}</text>'
        )

    svg.append('</svg>')
    svg = "".join(svg)
    svg = svg.replace("\n", "")

    return svg


def build_blast_reference_data(blast_file):
    """
    Return per-reference data for embedding in reports (no file writing)
    """
    ref_lengths, all_contig_aligns = parse_blast_for_ruler_continuous(blast_file)
    rows = extract_top_hit_table_with_alignment(blast_file)

    # --- group by reference ---
    ref_to_contigs = {}
    ref_to_table_rows = {}

    for block in all_contig_aligns:
        ref = block['ref']
        ref_clean = ref.split('.')[0]
        ref_to_contigs.setdefault(ref_clean, []).append(block)

    for row in rows:
        ref = row[1]
        ref_clean = ref.split('.')[0]
        ref_to_table_rows.setdefault(ref_clean, []).append(row)

    # --- build output ---
    result = {}

    for ref_id in ref_to_contigs:

        # 1. SVG (inline string, NOT file)
        svg_string = generate_reference_svg(
            ref_lengths[ref_id],
            ref_to_contigs[ref_id]
        )

        # 2. Alignment HTML per contig
        alignments_dict = {}
        for row in ref_to_table_rows.get(ref_id, []):
            query = row[0].split()[0].strip()
            alignment_block = row[-1]

            alignments_dict[query] = highlight_alignment_block(alignment_block)

        result[ref_id] = {
            "svg": svg_string,
            "alignments": alignments_dict,
            "table_html": output_alignments_html(ref_to_table_rows.get(ref_id, []))
        }


    return result
def wrap_sequence(seq, width=50):
    return "\n".join(seq[i:i+width] for i in range(0, len(seq), width))
    
def build_contig_orf_data(orf_fasta, contig_lengths=None):
    contigs = {}

    current_header = None
    current_seq = []

    with open(orf_fasta) as f:
        for line in f:
            line = line.strip()

            # --- New header ---
            if line.startswith(">"):
                # Save previous ORF before starting new one
                if current_header:
                    _store_orf(current_header, "".join(current_seq), contigs, contig_lengths)

                current_header = line[1:]
                current_seq = []

            else:
                current_seq.append(line)

        # Save last ORF
        if current_header:
            _store_orf(current_header, "".join(current_seq), contigs, contig_lengths)

    return contigs


def _store_orf(header, sequence, contigs, contig_lengths):
    # --- Extract contig ---
    contig = header.split("_ORF")[0]

    # --- Extract ORF ID ---
    orf_match = re.search(r"_ORF\.(\d+)", header)
    orf_id = orf_match.group(1) if orf_match else "?"

    # --- Extract coordinates ---
    coord_match = re.search(r"\[(\d+)-(\d+)\]\(([+-])\)", header)
    if not coord_match:
        return

    start, end, strand = coord_match.groups()

    # --- init contig ---
    if contig not in contigs:
        contigs[contig] = {"orfs": []}
        if contig_lengths:
            contigs[contig]["length"] = contig_lengths.get(contig)

    # --- add ORF ---
    contigs[contig]["orfs"].append({
        "id": orf_id,
        "start": int(start),
        "end": int(end),
        "strand": strand,
        "sequence": sequence,
#        "sequence_wrapped": wrap_sequence(sequence)
    })

def parse_contig_lengths(stats_file):
    df = pd.read_csv(stats_file, sep="\t")
    df.columns = [c.strip() for c in df.columns]

    # --- normalize contig name ---
    df["contig"] = df["qseqid"].str.split().str[0]

    return dict(zip(df["contig"], df["qlen"]))



def parse_hmmscan_domains(hmmscan_file):
    domains = {}

    with open(hmmscan_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue

            parts = line.split()

            target = parts[0]
            accession = parts[1]
            query = parts[3]

            # extract contig + ORF
            m = re.match(r"(.*)_ORF\.(\d+)", query)
            if not m:
                continue

            contig, orf_id = m.groups()

            # domain coordinates (ali coord)
            start = int(parts[17])
            end   = int(parts[18])

            desc = " ".join(parts[22:])
            i_evalue = parts[12]
            if float(i_evalue) > 1e-5:
                continue

            domains.setdefault(contig, {}).setdefault(orf_id, []).append({
                "name": target,
                "accession": accession,
                "start": start,
                "end": end,
                "desc": desc
            })

    return domains