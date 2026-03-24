#!/usr/bin/env python

import re
import os
import matplotlib.pyplot as plt
import html

# --- Utility: BLAST Parsing to get all reference & per-contig alignment boxes (one per contig/ref) ---
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
                    ref_lengths[current_ref] = int(m_len.group(1))
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

# --- Utility: Draw SVG ruler for reference+contigs ---
def plot_reference_ruler_from_blast(ref_lengths, contig_alignments, ref_id=None, svg_out="reference_contigs.svg"):
    for r_id, r_len in ref_lengths.items():
        if ref_id and r_id != ref_id: continue
        relevant_ctgs = [c for c in contig_alignments if c['ref'] == r_id]
        fig, ax = plt.subplots(figsize=(8, max(2, 0.25*(len(relevant_ctgs)+3))))
        y0 = 0
        ax.plot([0, r_len], [y0, y0], lw=10, color='#1267b2')
        for tick in range(0, r_len+1, 1000):
            ax.text(tick, y0-0.3, f"{tick//1000}k" if tick else "0k", ha="center", va="top", fontsize=10)
            ax.plot([tick, tick], [y0+0.2, y0-0.2], lw=1, color='black')
        spacing = 0.6
        for idx, ca in enumerate(relevant_ctgs):
            y = y0 + 0.8 + idx*spacing
            rect = plt.Rectangle((ca['ref_start'], y), ca['ref_end']-ca['ref_start'], 0.30, color='#e33')
            ax.add_patch(rect)
            ax.text((ca['ref_start']+ca['ref_end'])/2, y+0.19, ca["contig"], ha="center", va="bottom", fontsize=9)
        ax.set_xlim(-0.01*r_len, r_len*1.01)
        ax.set_ylim(y0-0.5, y0+1.2+len(relevant_ctgs)*spacing)
        ax.axis("off")
        fig.tight_layout()
        fig.savefig(svg_out, format="svg", bbox_inches="tight")
        plt.close(fig)

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
  <td style="color:#178cd2;"><a name="{query}">{query}</a></td>
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

# --- Index Page Generator ---
def write_reference_index_html(ref_ids, output_dir):
    html = ['<h2>References</h2>', '<ul>']
    for ref in ref_ids:
        html.append(f'<li><a href="{ref}.html">{ref}</a></li>')
    html.append('</ul>')
    with open(f"{output_dir}/index.html", "w") as fout:
        fout.write('\n'.join(html))

# --- Per-reference detail HTML page generator ---
def write_reference_html(ref_id, ref_len, contig_blocks, table_rows, output_dir):
    svg_filename = f"{output_dir}/{ref_id}_ruler.svg"
    plot_reference_ruler_from_blast({ref_id: ref_len}, contig_blocks, ref_id=ref_id, svg_out=svg_filename)
    html = f"""<html>
<head><title>Alignments for {ref_id}</title></head>
<body>
<div><a href="index.html">&larr; Back to References</a></div>
<h2>Reference: {ref_id}</h2>
<div style="text-align: center; margin-bottom:18px;">
  <img src="{ref_id}_ruler.svg" alt="Reference {ref_id} ruler" style="max-width:98%;">
</div>
"""
    html += output_alignments_html(table_rows, svg_filename=None)  # img already shown
    html += "</body></html>"
    with open(f"{output_dir}/{ref_id}.html", "w") as fout:
        fout.write(html)

def main():
    import argparse
    parser = argparse.ArgumentParser(description='Generate Per-Reference HTML reports from BLAST output.')
    parser.add_argument('-i', '--input', type=str, required=True, help='Input BLAST file')
    parser.add_argument('-o', '--outdir', type=str, required=True, help='Output HTML directory')
    args = parser.parse_args()

    output_dir = args.outdir
    os.makedirs(output_dir, exist_ok=True)

    ref_lengths, all_contig_aligns = parse_blast_for_ruler_continuous(args.input)
    rows = extract_top_hit_table_with_alignment(args.input)

    # group contig_blocks and table rows for each reference
    ref_to_contigs = {}
    ref_to_table_rows = {}
    for block in all_contig_aligns:
        ref_to_contigs.setdefault(block['ref'], []).append(block)
    for row in rows:
        ref = row[1]
        ref_to_table_rows.setdefault(ref, []).append(row)

    all_refs = sorted(ref_to_contigs.keys())
    write_reference_index_html(all_refs, output_dir)
    for ref in all_refs:
        write_reference_html(
            ref_id=ref,
            ref_len=ref_lengths[ref],
            contig_blocks=ref_to_contigs[ref],
            table_rows=ref_to_table_rows.get(ref, []),
            output_dir=output_dir,
        )

if __name__ == "__main__":
    main()