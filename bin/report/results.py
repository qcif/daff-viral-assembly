"""Define specific results from the analysis."""

import csv
from typing import Optional, Union, get_args, get_origin

from Bio import SeqIO

from .config import Config

config = Config()

COLLECT_RANKS = [
    'acellular root',
    'domain',
    'kingdom',
    'phylum',
    'class',
    'order',
    'family',
    'genus',
    'species',
]


class FLAGS:
    SUCCESS = 'success'
    WARNING = 'warning'
    DANGER = 'danger'
    NONE = 'secondary'


def colour_to_bs_class(colour: str) -> str:
    """Convert a colour string to a Bootstrap class."""
    mapping = {
        'green': FLAGS.SUCCESS,
        'yellow': FLAGS.WARNING,
        'red': FLAGS.DANGER,
    }
    return mapping.get(colour.lower(), FLAGS.NONE)


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
        ('fastq_1', str),
        ('fastq_2', str),
        ('sample_information', str),
        ('sample_type', str),
        ('sample_receipt_date', str),
        ('storage_location', str),
    ]


class RunQC(AbstractDataRow):
    """Report the quality control outcomes."""

    COLUMNS = [
        ('raw_reads', int),
        ('quality_filtered_reads', int),
        ('mean_raw_read_length', int),
        ('mean_filtered_read_length', int),
        ('gc_content', float),
        ('rRNA_cleaned_reads', int),
        ('phix_cleaned_reads', int),
        ('percent_qfiltered', float),
        ('percent_cleaned', float),
        ('raw_reads_flag', Optional[str]),
        ('qfiltered_reads_flag', Optional[str]),
        ('QC_FLAG', str),
    ]

    @property
    def flag(self):
        return colour_to_bs_class(self.QC_FLAG)

        #! Check we don't want a calculated flag?
        # raw_threshold = self.raw_reads > config.CRITERIA.MIN_RAW_READS
        # qfiltered_threshold = (
        #     (self.percent_cleaned or 0 / 100 * self.raw_reads)
        #     > config.CRITERIA.MIN_FILTERED_READS)
        # if qfiltered_threshold:
        #     if raw_threshold:
        #         return FLAGS.SUCCESS
        #     return FLAGS.WARNING
        # return FLAGS.DANGER

    @property
    def html_file(self):
        return config.run_qc_html_file


class AbstractResultRows:
    """A result composed of a series of rows with defined columns names.

    This allows for easy parsing and rendering of different tabular data by
    setting the COLUMNS and COLUMN_METADATA attributes, which should describe
    the table schema.

    COLUMNS: list of column names in order
    COLUMN_METADATA: dict of column name to metadata dict, which should include
        'type': str indicating type for casting (e.g., 'int', 'float', ...) and
        and other arbitrary metadata that might be used in the subclass.

    """
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

    @classmethod
    def from_csv(cls, path, delimiter='\t'):
        """Return the blast hits as a dict."""
        with path.open() as f:
            reader = csv.DictReader(f, delimiter=delimiter)
            return cls(reader)

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
        try:
            if not type_str:
                return value
            if type_str == 'int':
                num = int(float(value))
                return f"{num:,}"
            if type_str == 'float':
                if value.lower() in ('na', 'nan', 'n/a'):
                    return None
                return float(value)
            if type_str == 'scientific' and 'e' in value:
                num = float(value)
                return f"{num:.2e}"
            if type_str == 'bool':
                return value.lower() in ['true', '1', 'yes', 'y']
        except ValueError:
            return None
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
        self.columns_secondary_display = [
            c for c in self.columns_display
            if self.COLUMN_METADATA[c]['secondary_display']
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
            # TODO: define this
            return 'success'

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
                for colname in self.COLUMNS[3:]:
                    row[colname] = '-'


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


class AbstractTaxonomicClassificationResult(AbstractResultRows):
    """Base class for taxonomic classification results (Kraken, Kaiju)."""

    def __init__(self, *args):
        super().__init__(*args)
        self.set_taxonomy_columns()
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
            # TODO: define this
            return 'success'

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
            if not row['taxon_id'] or row['taxon_id'] == '0':
                for colname in self.COLUMNS[3:]:
                    row[colname] = '-'

    def set_taxonomy_columns(self):
        """Set taxonomy columns based on full_lineage field."""
        for row in self.rows:
            lineage = row.get('full_lineage', '').split(';')
            ranks = [
                rank.lower()
                for rank in row.get('full_lineage_ranks', '').split(';')
            ]
            taxonomy = dict(zip(ranks, lineage))
            for rank in COLLECT_RANKS:
                row[rank] = taxonomy.get(rank, '-')


class KrakenResults(AbstractTaxonomicClassificationResult):
    COLUMN_METADATA = _csv_to_dict(config.SCHEMA.KRAKEN_FIELD_CSV)
    COLUMNS = list(COLUMN_METADATA.keys())


class KaijuResults(AbstractTaxonomicClassificationResult):
    COLUMN_METADATA = _csv_to_dict(config.SCHEMA.KAIJU_FIELD_CSV)
    COLUMNS = list(COLUMN_METADATA.keys())


class MappingResults(AbstractResultRows):
    COLUMN_METADATA = _csv_to_dict(config.SCHEMA.MAPPING_FIELD_CSV)
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
            # TODO: define this
            return 'success'

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
                for colname in self.COLUMNS[3:]:
                    row[colname] = '-'


class SummaryResults(AbstractResultRows):
    COLUMN_METADATA = _csv_to_dict(config.SCHEMA.SUMMARY_FIELD_CSV)
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


class NovelVirusResults(AbstractResultRows):
    COLUMN_METADATA = _csv_to_dict(config.SCHEMA.NOVEL_VIRUS_FIELD_CSV)
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
        for row in self.rows:
            if row.get('taxonomy'):
                row['taxonomy'] = [
                    taxon.strip()
                    for taxon in row.get('taxonomy', '').split(';')
                    if taxon.strip()
                ]

    def taxonomy_html_for_row(self, row_ix):
        """Return the taxonomy as multiline HTML string."""
        row = self.rows[row_ix]
        return '<br>'.join(row.get('taxonomy', []))
