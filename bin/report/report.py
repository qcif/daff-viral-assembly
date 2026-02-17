"""Render a workflow report from output files."""

import csv
import json
import logging
import os
from datetime import datetime
from pathlib import Path

import yaml
from jinja2 import Environment, FileSystemLoader

from . import config
# from .bam import render_bam_html
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

    # if len(context['blast_hits']):
    #     render_bam_html()


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
        'flags': config.flags,
        'blast_passed': config.blast_passed,
        'filter_keywords': config.filter_keywords.read_text().splitlines(),
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
    return {
        'count': len(blast_hits.positive_hits),
        'percent': 'NA' if not contigs_fasta else round(
            100 * len(blast_hits.positive_hits) / len(contigs_fasta)
        ),
    }
