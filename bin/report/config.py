"""Configuration of the report module."""

import csv
import os
from datetime import datetime
from functools import cached_property
from pathlib import Path

import yaml

from .utils import path_safe

PARENT_DIR = Path(__file__).parent
ROOT_DIR = Path(__file__).parents[2].resolve()
REPO_URL = 'https://github.com/qcif/daff-viral-assembly'
WF_NAME = 'ViroGen'


class Config:

    TIMESTAMP_FILE = '*_start_timestamp.txt'
    VERSIONS_PATH = ROOT_DIR / 'versions.yml'
    DEFAULT_PARAMS_PATH = ROOT_DIR / 'params/default_params.yml'
    FLAGS_CSV = ROOT_DIR / 'flags.csv'

    class SCHEMA:
        BLAST_HITS_FIELD_CSV = PARENT_DIR / 'schema/blast_fields.csv'
        KRAKEN_FIELD_CSV = PARENT_DIR / 'schema/kraken_fields.csv'
        KRAKEN_KAIJU_FIELD_CSV = PARENT_DIR / 'schema/kraken_kaiju_fields.csv'
        MAPPING_FIELD_CSV = PARENT_DIR / 'schema/mapping_fields.csv'
        SUMMARY_FIELD_CSV = PARENT_DIR / 'schema/summary_fields.csv'
        NOVEL_VIRUS_FIELD_CSV = PARENT_DIR / 'schema/novel_virus_fields.csv'

    class OUTPUTS:
        REPORT_FILE_TEMPLATE = '{sample_id}_report.html'
        BAM_HTML_FILE_TEMPLATE = '{sample_id}_bam-alignment.html'

    class REPORT:
        TITLE = "Viral genome assembly report"
        SUBTITLE = (
            "Results generated through the"
            f' <a href="{REPO_URL}" target="_blank">'
            f'{WF_NAME} NextFlow workflow</a>.')

    class CRITERIA:
        MIN_RAW_READS = 8000000  # ! confirm with DAFF
        MIN_FILTERED_READS = 2500000  # ! confirm with DAFF

    @property
    def default_params(self) -> dict[str, str]:
        """Return dict of default workflow parameters."""
        with self.DEFAULT_PARAMS_PATH.open() as f:
            return yaml.safe_load(f)

    @property
    def result_dir(self) -> Path:
        """Ensure that result dir is propagated between Config instances."""
        return Path(os.environ['RESULT_DIR'])

    @property
    def report_path(self) -> Path:
        path = self.result_dir / self.OUTPUTS.REPORT_FILE_TEMPLATE.format(
            sample_id=path_safe(self.sample_id),
        )
        return path.absolute()

    @property
    def bam_html_path(self) -> Path:
        path = self.result_dir / self.OUTPUTS.BAM_HTML_FILE_TEMPLATE.format(
            sample_id=path_safe(self.sample_id),
        )
        return path.absolute()

    @property
    def metadata_path(self) -> Path:
        return self.result_dir / self.METADATA_FILE

    @property
    def run_qc_path(self) -> Path:
        return self._get_file_by_pattern('run_qc_report_*.txt')

    @property
    def run_qc_html_file(self) -> str:
        return self._get_file_by_pattern('run_qc_report_*.html').name

    @property
    def bam_path(self) -> Path:
        return self._get_file_by_pattern("*.bam")

    @property
    def bai_path(self) -> Path:
        return self._get_file_by_pattern("*.bai")

    @property
    def contigs_fasta_path(self) -> Path:
        # ! confirm?
        return self._get_file_by_pattern("*_scaffolds.fasta")

    @property
    def reference_fasta_path(self) -> Path:
        return self._get_file_by_pattern("*_ref_sequences_clustered.fasta")

    @property
    def blast_hits_path(self) -> Path:
        return self._get_file_by_pattern(
            "*_megablast_top_viral_hits_with_contigs.txt")

    @property
    def kraken_hits_path(self) -> Path:
        return self._get_file_by_pattern("*_kraken_summary.txt")

    @property
    def kaiju_hits_path(self) -> Path:
        return self._get_file_by_pattern("*_kaiju_summary.txt")

    @property
    def ref_mapping_path(self) -> Path:
        return self._get_file_by_pattern(
            "*_reference_with_cov_stats_final.txt")

    @property
    def summary_path(self) -> Path:
        return self._get_file_by_pattern("*_summary_viral_results.tsv")

    @property
    def novel_viruses_path(self) -> Path:
        return self._get_file_by_pattern("*_novel_virus_candidates.tsv")

    # @property
    # def PLACEHOLDER(self) -> Path:
    #     return self._get_file_by_pattern("*.*")

    @property
    def blast_passed(self) -> bool:
        """Check if BLAST was successful."""
        try:
            path = self._get_file_by_pattern("*_blast_status.txt")
        except FileNotFoundError:
            return True  # If no file written, assume it passed
        if not path.exists():
            return True
        return 'fail' not in path.read_text().lower()

    @cached_property
    def sample_id(self) -> str:
        """Return the sample ID from the result directory."""
        bam_path = self._get_file_by_pattern('*.bam')
        return bam_path.name.split('_ref_aln.')[0]

    @cached_property
    def start_time(self) -> str:
        """Return the timestamp of the start of the workflow."""
        timestamp_path = self._get_file_by_pattern(
            self.TIMESTAMP_FILE,
            ignore_missing=True,
        )
        if timestamp_path and timestamp_path.exists():
            return datetime.strptime(
                timestamp_path.read_text().strip(),
                '%Y%m%d%H%M%S',
            )
        return None

    @cached_property
    def flags(self) -> list[dict[str, str]]:
        """Load flag definitions from a CSV file."""
        if not self.FLAGS_CSV.exists():
            return []
        with self.FLAGS_CSV.open() as f:
            reader = csv.DictReader(f)
            flags = list(reader)
        return flags

    def load(self, result_dir: Path):
        """Read the results from the result directory."""
        os.environ['RESULT_DIR'] = str(result_dir)

    def _get_file_by_pattern(
        self,
        file_pattern: str,
        ignore_missing=False,
    ) -> Path:
        paths = list(self.result_dir.glob(file_pattern, case_sensitive=False))
        if paths:
            return paths[0]
        if ignore_missing:
            return None
        raise FileNotFoundError(
            f'No file matching pattern: {self.result_dir / file_pattern}'
        )
