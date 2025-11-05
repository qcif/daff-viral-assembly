"""Define specific results from the analysis."""

import base64
import csv
from typing import Optional, Union, get_args, get_origin

from Bio import SeqIO

from .config import Config

config = Config()


def _csv_to_dict(csv_path, index_col='colname'):
    with open(csv_path) as f:
        reader = csv.reader(f)
        next(reader)  # skip header line
        ordered_colnames = [
            row[0].strip()
            for row in reader
        ]
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        data = {
            row[index_col].strip(): {
                colname: value
                for colname, value in row.items()
            }
            for row in reader
        }
        return {
            colname: data[colname]
            for colname in ordered_colnames
        }


class FLAGS:
    SUCCESS = 'success'
    WARNING = 'warning'
    DANGER = 'danger'
    NONE = 'secondary'


class AbstractDataRow:

    COLUMNS = []

    def __init__(self, row):
        for colname, _type in self.COLUMNS:
            raw_value = row.get(colname, None)
            if raw_value is None:
                value = None
            else:
                # If it's an Optional (i.e., Union[..., NoneType])
                origin = get_origin(_type)
                if origin is Union:
                    # Extract actual type from Optional[Type]
                    allowed_types = [
                        t for t in get_args(_type)
                        if t is not type(None)
                    ]
                    if allowed_types:
                        value = allowed_types[0](raw_value.strip())
                    else:
                        value = None
                else:
                    value = _type(raw_value.strip())

            setattr(self, colname, value)

    def to_json(self):
        return {
            colname: getattr(self, colname)
            for colname, _ in self.COLUMNS
        }


class Metadata(AbstractDataRow):
    """Report the sample metadata."""

    COLUMNS = [
        ('fastq_path', str),
        ('target_size', int),
        ('target_organism', str),
        ('target_gene', str),
        ('fwd_primer', Optional[str]),
        ('rev_primer', Optional[str]),
        ('test', Optional[str]),
        ('method', Optional[str]),
    ]


class RunQC(AbstractDataRow):
    """Report the quality control outcomes."""

    COLUMNS = [
        ('raw_reads', int),
        ('processed_reads', int),
        ('percent_processed', float),
        ('raw_reads_flag', str),
        ('processed_flag', str),
    ]

    @property
    def flag(self):
        raw_threshold = self.raw_reads > config.CRITERIA.MIN_RAW_READS
        qfiltered_threshold = (
            (self.processed_reads or 0)
            > config.CRITERIA.MIN_FILTERED_READS)
        if qfiltered_threshold:
            if raw_threshold:
                return FLAGS.SUCCESS
            return FLAGS.WARNING
        return FLAGS.DANGER

    @property
    def nanoplot_raw_html_base64(self):
        content_bytes = config.nanoplot_raw_html_path.read_bytes()
        return base64.b64encode(content_bytes).decode()

    @property
    def nanoplot_filtered_html_base64(self):
        content_bytes = config.nanoplot_filtered_html_path.read_bytes()
        return base64.b64encode(content_bytes).decode()

    @property
    def html_file(self):
        return config.run_qc_html_file


class AbstractResultRows:
    """A result composed of a series of rows with defined columns names."""
    COLUMNS = []
    COLUMN_METADATA = {}

    def __init__(self, rows):
        self.rows = self._parse_rows(rows)

    def __len__(self):
        return len(self.rows)

    def __iter__(self):
        return iter(self.rows)

    def __getitem__(self, index):
        return self.rows[index]

    def _parse_rows(self, rows):
        return [
            {
                colname: (
                    self._cast(
                        row[colname].strip(),
                        self.COLUMN_METADATA.get(colname, {}).get('type'),
                    )
                )
                for colname in self.COLUMNS
                if colname in row
            }
            for row in rows
        ]

    def _cast(self, value, type_str):
        if not type_str:
            return value
        if type_str == 'int':
            num = int(float(value))
            return f"{num:,}"
        if type_str == 'float':
            return float(value)
        if type_str == 'scientific' and 'e' in value:
            num = float(value)
            return f"{num:.2e}"
        return value

    def to_json(self):
        return self.rows


class BlastHits(AbstractResultRows):
    COLUMN_METADATA = _csv_to_dict(config.SCHEMA.BLAST_HITS_FIELD_CSV)
    COLUMNS = list(COLUMN_METADATA.keys())

    def __init__(self, *args):
        super().__init__(*args)
        self.columns_display = [
            c for c in self.COLUMN_METADATA
            if self.COLUMN_METADATA[c]['label']
        ]
        self.columns_primary_display = [
            c for c in self.columns_display
            if self.COLUMN_METADATA[c]['primary_display']
        ]
        self.set_bs_class()
        self.set_null_rows()

    @property
    def positive_hits(self):
        return [
            row for row in self.rows
            if row.get('sacc') not in [None, '', '0', '-']
        ]

    def set_bs_class(self):
        def _get_bs_class(row):
            # Extract numeric scores (expected in the range 0 to 1)
            conf_score = row.get('NORMALISED_CONF_SCORE', 1)
            sacc = row.get('sacc')
 
            # Define thresholds (customize as needed)
            if conf_score == 0 and sacc in [None, '', '0', '-']:
                return 'secondary'   # grey
            elif conf_score < 0.5:
                return 'danger'      # red
            elif conf_score < 0.8:
                return 'warning'     # orange
            else:
                return 'success'     # green

        self.rows = [
            {
                **row,
                'bs_class': _get_bs_class(row),
            }
            for row in self.rows
        ]

    def set_null_rows(self):
        """Set rows with no hits to have a null value."""
        for row in self.rows:
            if not row['sacc'] or row['sacc'] == '0':
                for colname in self.COLUMNS[3:41] + self.COLUMNS[49:50]:
                    row[colname] = '-'


class BlastHitsPolished(AbstractResultRows):
    COLUMNS = [
        'qseqid',
        'sgi',
        'sacc',
        'length',
        'nident',
        'pident',
        'mismatch',
        'gaps',
        'gapopen',
        'qstart',
        'qlen',
        'qend',
        'sstart',
        'send',
        'slen',
        'sstrand',
        'evalue',
        'bitscore',
        'qcovhsp',
        'stitle',
        'staxids',
        'qseq',
        'sseq',
        'sseqid',
        'qcovs',
        'qframe',
        'sframe',
        'species',
        'broad_taxonomic_category',
    ]


class ConsensusFASTA:

    def __init__(self, fasta_path):
        self.fasta_path = fasta_path
        self.records = list(SeqIO.parse(fasta_path, 'fasta'))

    def __len__(self):
        return len(self.records)

    def __bool__(self):
        return bool(len(self.records))

    def __iter__(self):
        return iter(self.records)

    def __str__(self):
        return '\n\n'.join(
            f">{seq.id}\n"
            + self._wrap(seq.seq)
            for seq in self.records
        )

    def _wrap(self, seq, width=80):
        seq_str = str(seq)
        return '\n'.join(
            seq_str[i:i + width]
            for i in range(0, len(seq), width)
        )

    def to_json(self):
        return {
            seq.id: str(seq.seq)
            for seq in self.records
        }
