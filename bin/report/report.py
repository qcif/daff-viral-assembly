"""Render a workflow report from output files."""

import csv
import json
import logging
import os
import re
import sys
from datetime import datetime
from pathlib import Path

import yaml
from jinja2 import Environment, FileSystemLoader
#from .blast_viz import build_blast_reference_data, parse_contig_lengths, build_contig_orf_data, _store_orf
from . import blast_viz
from . import config
from .bam import render_bam_html
from .results import (
    BlastHits,
    ConsensusFASTA,
    KaijuResults,
    KrakenResults,
    MappingResults,
    Metadata,
    RunQC,
    SummaryResults,
    NovelVirusResults,
)
from .filters.css_hash import css_hash
from .utils import get_img_src, serialize
from .taxa import parse_kraken_taxonomy

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)
config = config.Config()

TEMPLATE_DIR = Path(__file__).parent / 'templates'
STATIC_DIR = Path(__file__).parent / 'static'
EXCLUDE_JS = [
    'igv-',
]


def render(
    result_dir: Path,
    samplesheet_file: Path,
    default_params_file: Path,
    params_file: Path,
    versions: Path,
    analyst_name: str = None,
    facility: str = None
):
    """Render to HTML report to the configured output directory."""
    config.load(result_dir)
    j2 = Environment(loader=FileSystemLoader(TEMPLATE_DIR))
    j2.filters['css_hash'] = css_hash
    j2.filters['html_id'] = lambda s: re.sub(r'[^A-Za-z0-9_-]', '_', str(s))
    template = j2.get_template('index.html')
    context = _get_report_context(
        samplesheet_file,
        default_params_file,
        params_file,
        versions,
        analyst_name,
        facility,
    )
    path = config.result_dir / 'report_context.json'
    with path.open('w') as f:
        logger.info(f"Writing report context to {path}")
        json.dump(context, f, indent=2, default=serialize)

    static_files = _get_static_file_contents()
    rendered_html = template.render(**context, **static_files)

    path = config.report_path
    with open(path, 'w') as f:
        f.write(rendered_html)
    logger.info(f"HTML document written to {path}")

    if len(context['blast_hits']) and not os.getenv('SKIP_BAM_RENDER'):
        render_bam_html()


def _get_static_file_contents():
    """Return the static files content as strings."""
    static_files = {}
    for root, _, files in os.walk(STATIC_DIR):
        root = Path(root)
        if root.name == 'css':
            static_files['css'] = [
                f'/* {f} */\n' + (root / f).read_text()
                for f in sorted(files)
            ]
        elif root.name == 'js':
            static_files['js'] = [
                f'/* {f} */\n' + (root / f).read_text()
                for f in sorted(files)
                if not any(
                    phrase in f
                    for phrase in EXCLUDE_JS
                )
            ]
        elif root.name == 'img':
            static_files['img'] = {
                f: get_img_src(root / f)
                for f in sorted(files)
            }

    return {'static': static_files}


def _get_report_context(
    samplesheet_file: Path,
    default_params_file: Path,
    params_file: Path,
    versions: Path,
    analyst_name: str,
    facility: str,
) -> dict:
    """Build the context for the report template."""
    blast_hits = BlastHits.from_csv(config.blast_hits_path)
    contigs_fasta = ConsensusFASTA(config.contigs_fasta_path)
    blast_reference_data = {}
    if config.blast_output_path and config.blast_output_path.exists():
        logger.info(f"Building BLAST reference data from {config.blast_output_path}")
        blast_reference_data = blast_viz.build_blast_reference_data(config.blast_output_path)
    contig_lengths = blast_viz.parse_contig_lengths(config.blast_hits_path)
    orf_data = blast_viz.build_contig_orf_data(config.orf_fasta_path, contig_lengths)
    
    domain_data = blast_viz.parse_hmmscan_domains(config.hmmscan_output_path)

    for contig_id, data in orf_data.items():

        orfs = data.get("orfs", [])
        if not orfs:
            continue

        length = data.get("length")
        domains_for_contig = domain_data.get(contig_id, {})

        for orf in orfs:
            orf["domains"] = domains_for_contig.get(orf["id"], [])

        if length:
            data["svg"] = generate_contig_orf_svg(contig_id, length, orfs)

    return {
        'title': config.REPORT.TITLE,
        'subtitle_html': config.REPORT.SUBTITLE,
        'sample_id': config.sample_id,
        'analyst_name': analyst_name or '-',
        'facility': facility or '-',
        'versions': _get_versions(versions),
        'start_time': _get_start_time(),
        'end_time': datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        'wall_time': _get_walltime(),
        'metadata': _get_metadata(samplesheet_file),
        'parameters': _get_parameters(default_params_file, params_file),
        'run_qc': _get_run_qc(),
        'run_qc_raw_fastqc_filename': config.run_qc_raw_fastqc_path.name,
        'run_qc_clean_fastqc_filename': config.run_qc_clean_fastqc_path.name,
        'read_distribution_chart_src': get_img_src(
            config.read_distribution_chart_path
        ) if config.read_distribution_chart_path.exists() else None,
        'bam_html_file': config.bam_html_path.name,
        'contigs_fasta': contigs_fasta,
        'blast_hits': blast_hits,
        'blast_stats': _calculate_blast_stats(
            blast_hits,
            contigs_fasta,
        ),
        'kraken': KrakenResults.from_csv(config.kraken_hits_path),
        'kaiju': KaijuResults.from_csv(config.kaiju_hits_path),
        'mapping': MappingResults.from_csv(config.ref_mapping_path),
        'summary': SummaryResults.from_csv(config.summary_path),
        'novel_viruses': NovelVirusResults.from_csv(config.novel_viruses_path),
        'kraken_taxa_by_kingdom': parse_kraken_taxonomy(
            config.kraken_hits_path,
            group_kingdoms=True,
        ),
        'kaiju_taxa_by_kingdom': parse_kraken_taxonomy(
            config.kaiju_hits_path,
            group_kingdoms=True,
        ),
        'flags': config.flags,
        'blast_passed': config.blast_passed,
        'filter_keywords': config.filter_keywords.read_text().splitlines(),
        'blast_reference_data': blast_reference_data,
        'contig_orf_data': orf_data,
    }


def _get_start_time():
    if not config.start_time:
        return None
    return config.start_time.strftime("%Y-%m-%d %H:%M:%S")


def _get_walltime():
    """Return wall time since start of the workflow.
    Returns a dict of hours, minutes, seconds.
    """
    if not config.start_time:
        return None
    seconds = (datetime.now() - config.start_time).total_seconds()
    hours, remainder = divmod(seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    return {
        'hours': str(int(hours)).zfill(2),
        'minutes': str(int(minutes)).zfill(2),
        'seconds': str(int(seconds)).zfill(2),
    }


def _get_metadata(samplesheet_file: Path):
    """Return the metadata as a dict."""
    with open(samplesheet_file) as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row['sample'] == config.sample_id:
                return Metadata(row)


def _get_default_parameters(default_params_file: Path) -> dict[str, str]:
    """Return the default_parameters as a dict."""
    if not default_params_file or not default_params_file.exists():
        return {}
    with default_params_file.open() as f:
        return yaml.safe_load(f)


def _get_parameters(
    default_params_file: Path,
    params_file: Path,
) -> dict[str, dict[str, str]]:
    """Return the parameters as a dict."""
    if not default_params_file or not default_params_file.exists():
        return {}
    defaults = _get_default_parameters(default_params_file)
    with params_file.open() as f:
        user_params = yaml.safe_load(f)
    return {
        k: {
            'default': defaults[k],
            'user': user_params.get(k, None),
        } for k in defaults
    }


def _get_versions(versions: Path) -> dict[str, str]:
    """Return dict of program versions used in the workflow."""
    if not versions or not versions.exists():
        return {}
    with versions.open() as f:
        return yaml.safe_load(f)


def _get_run_qc() -> dict:
    """Return the runs stats as a dict.

    Columns:
    - Sample
    - raw_reads
    - quality_filtered_reads
    - mean_raw_read_length
    - mean_filtered_read_length
    - gc_content
    - rRNA_cleaned_reads
    - phix_cleaned_reads
    - percent_qfiltered
    - percent_cleaned
    - raw_reads_flag
    - qfiltered_reads_flag
    - QC_FLAG

    """
    with config.run_qc_path.open() as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if row['Sample'] == config.sample_id:
                return RunQC(row)
    return {}


def _calculate_blast_stats(
    blast_hits: BlastHits,
    contigs_fasta: ConsensusFASTA,
) -> dict:
    percent_hits = (
        100 * len(blast_hits.positive_hits)
        / len(contigs_fasta)
    )
    return {
        'count': len(blast_hits.positive_hits),
        'percent': 'NA' if not contigs_fasta else float(
            f'{percent_hits:.2g}'
        ),
    }

def generate_contig_orf_svg(contig_id, contig_length, orfs):
    import re

    width = 800
    scale = width / contig_length

    orfs = assign_tracks(orfs)

    track_height = 35
    contig_y = 30
    orf_gap = 20   # space between contig and first ORF

    y_offset = contig_y + orf_gap

    domain_height = 6
    domain_offset = 16  # below ORF arrow

    safe_contig_id = re.sub(r'[^A-Za-z0-9_-]', '_', contig_id)

    # dynamic height
    n_tracks = max([orf['track'] for orf in orfs], default=0) + 1
    height = y_offset + n_tracks * track_height + 20

    svg = [f'<svg width="{width}" height="{height}" xmlns="http://www.w3.org/2000/svg">']

    # --- contig line ---
    svg.append(
        f'<line x1="0" y1="{contig_y}" x2="{width}" y2="{contig_y}" stroke="#1267b2" stroke-width="8"/>'
    )

    # --- ORFs ---
    for orf in orfs:
        x1 = orf['start'] * scale
        x2 = orf['end'] * scale
        y = y_offset + orf['track'] * track_height

        color = "#2ca02c" if orf['strand'] == '+' else "#d62728"
        arrow_size = 6

        if orf['strand'] == '+':
            points = [
                (x1, y),
                (x2 - arrow_size, y),
                (x2, y + 5),
                (x2 - arrow_size, y + 10),
                (x1, y + 10)
            ]
        else:
            points = [
                (x2, y),
                (x1 + arrow_size, y),
                (x1, y + 5),
                (x1 + arrow_size, y + 10),
                (x2, y + 10)
            ]

        # --- clickable ORF ---
        svg.append(
            f'''
        <a xlink:href="javascript:void(0)"
        data-bs-toggle="modal"
        data-bs-target="#orfModal_{safe_contig_id}_{orf['id']}">
        <polygon points="{" ".join(f"{px},{py}" for px, py in points)}"
                fill="{color}" />
        </a>
        '''
        )

        # --- ORF label ---
        svg.append(
            f'<text x="{(x1+x2)/2}" y="{y-2}" font-size="9" text-anchor="middle">{orf["id"]}</text>'
        )

        # --- domains ---
        orf_length = max(1, orf["end"] - orf["start"] + 1)  # avoid /0

        domains = assign_domain_tracks(orf.get("domains") or [])
        #for domain in orf.get("domains", []):
        for domain in domains:
            d1 = x1 + ((domain["start"])     / orf_length) * (x2 - x1)
            d2 = x1 + ((domain["end"])       / orf_length) * (x2 - x1)

            y_domain = y + domain_offset + domain["track"] * (domain_height + 2)

            svg.append(
                f'''
            <g>
            <title>{domain["name"]} ({domain["accession"]})</title>
            <rect x="{d1}" y="{y_domain}"
                    width="{d2 - d1}" height="{domain_height}"
                    fill="#9467bd" />
            </g>
            '''
            )

    svg.append('</svg>')
    return "".join(svg)

def assign_tracks(orfs):
    tracks = []

    for orf in sorted(orfs, key=lambda x: x['start']):
        placed = False

        for i, track in enumerate(tracks):
            # check overlap with last ORF in track
            if orf['start'] > track[-1]['end']:
                track.append(orf)
                orf['track'] = i
                placed = True
                break

        if not placed:
            orf['track'] = len(tracks)
            tracks.append([orf])

    return orfs


def assign_domain_tracks(domains):
    domains = sorted(domains, key=lambda d: d["start"])
    tracks = []

    for d in domains:
        placed = False

        for i, track in enumerate(tracks):
            if d["start"] > track[-1]["end"]:
                track.append(d)
                d["track"] = i
                placed = True
                break

        if not placed:
            tracks.append([d])
            d["track"] = len(tracks) - 1

    return domains