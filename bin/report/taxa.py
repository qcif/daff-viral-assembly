"""Parse taxonomic classifications from Kraken summary files."""

import csv
from collections import defaultdict
from pathlib import Path
from typing import Union


TOP_N_TAXA = 10
TAXONOMIC_RANKS = [
    'domain',
    'kingdom',
    'phylum',
    'class',
    'order',
    'family',
    'genus',
    'species',
]
KINGDOM_GROUPS = [
    'plant',
    'animal',
    'bacteria',
    'fungi',
    'virus',
    'other',
]
PLANT_MARKERS = ['viridiplantae', 'plantae']
FUNGI_MARKERS = ['fungi', 'ascomycota', 'basidiomycota']
ANIMAL_MARKERS = ['metazoa', 'chordata', 'vertebrata']


def parse_kraken_taxonomy(
    kraken_summary_path: Union[str, Path],
    group_kingdoms: bool = False,
) -> dict[str, dict[str, dict[str, int]]]:
    """Parse Kraken summary file and enumerate taxa by rank.

    Args:
        kraken_summary_path: Path to the Kraken summary TSV file
        group_kingdoms: If True, group species by kingdom categories
                       (plant, animal, bacteria, fungi, virus, other)
                       instead of returning all taxonomic ranks

    Returns:
        When group_kingdoms=False (default):
            A nested dictionary structure:
            {
                rank: {
                    taxon: {
                        'taxon_count': <int>,  # number of rows
                        'read_count': <int>,   # sum of reads
                    }
                }
            }
            Taxa within each rank are ordered by descending read abundance.

        When group_kingdoms=True:
            A dictionary structure grouping taxa by kingdom:
            {
                kingdom_group: {
                    'read_count': <int>,   # total reads for kingdom
                    'taxon_count': <int>,  # total taxa for kingdom
                    'species': {
                        species_name: {
                            'taxon_count': <int>,  # number of rows
                            'read_count': <int>,   # sum of reads
                        }
                    },
                    'genus': {
                        genus_name: {
                            'taxon_count': <int>,  # number of rows
                            'read_count': <int>,   # sum of reads
                        }
                    },
                    'family': {
                        family_name: {
                            'taxon_count': <int>,  # number of rows
                            'read_count': <int>,   # sum of reads
                        }
                    }
                }
            }
            Taxa within each rank are limited to the top N (default 10) by
            read abundance, with remaining taxa grouped under 'other'.

    Example:
        >>> result = parse_kraken_taxonomy('kraken_summary.txt')
        >>> result['domain']['Eukaryota']['read_count']
        12345678
        >>> grouped = parse_kraken_taxonomy('kraken_summary.txt',
        ...                                  group_kingdoms=True)
        >>> grouped['plant']['read_count']
        52345678
        >>> grouped['plant']['species']['Citrus limon']['read_count']
        5180323
        >>> grouped['plant']['genus']['Citrus']['read_count']
        8234567
        >>> grouped['plant']['family']['Rutaceae']['read_count']
        12456789
    """
    kraken_summary_path = Path(kraken_summary_path)

    if group_kingdoms:
        return _parse_by_kingdom_groups(kraken_summary_path)
    else:
        return _parse_by_taxonomic_ranks(kraken_summary_path)


def _parse_by_taxonomic_ranks(
    kraken_summary_path: Path,
) -> dict[str, dict[str, dict[str, int]]]:
    """Parse Kraken summary by all taxonomic ranks.

    Args:
        kraken_summary_path: Path to the Kraken summary TSV file

    Returns:
        Dictionary organized by taxonomic ranks
    """
    # Initialize nested defaultdicts to accumulate data
    taxa_data = {rank: defaultdict(lambda: {'taxon_count': 0, 'read_count': 0})
                 for rank in TAXONOMIC_RANKS}

    with kraken_summary_path.open() as f:
        reader = csv.DictReader(f, delimiter='\t')

        for row in reader:
            full_lineage = row.get('full_lineage', '')
            full_lineage_ranks = row.get('full_lineage_ranks', '')
            reads = int(row.get('reads', 0))

            if not full_lineage or not full_lineage_ranks:
                continue

            # Split the lineage into taxonomic levels
            lineage_parts = [
                part.strip()
                for part in full_lineage.split(';')
                if part.strip()
            ]
            rank_parts = [
                part.strip()
                for part in full_lineage_ranks.split(';')
                if part.strip()
            ]

            # Map lineage parts to ranks using full_lineage_ranks
            rank_assignments = _assign_ranks_to_lineage(
                lineage_parts, rank_parts
            )

            # Accumulate counts for each rank
            for rank, taxon in rank_assignments.items():
                if rank in taxa_data:
                    taxa_data[rank][taxon]['taxon_count'] += 1
                    taxa_data[rank][taxon]['read_count'] += reads

    # Convert defaultdicts to regular dicts and sort by read_count (desc)
    result = {}
    for rank in TAXONOMIC_RANKS:
        sorted_taxa = dict(
            sorted(
                taxa_data[rank].items(),
                key=lambda x: x[1]['read_count'],
                reverse=True
            )
        )
        result[rank] = sorted_taxa

    return result


def _parse_by_kingdom_groups(
    kraken_summary_path: Path,
) -> dict[str, dict[str, dict[str, int]]]:
    """Parse Kraken summary grouping taxa by kingdom categories.

    Args:
        kraken_summary_path: Path to the Kraken summary TSV file

    Returns:
        Dictionary organized by kingdom groups with totals:
        {
            kingdom_group: {
                'read_count': <int>,     # total reads for kingdom
                'taxon_count': <int>,    # total taxa for kingdom
                'species': {
                    species_name: {
                        'taxon_count': <int>,
                        'read_count': <int>,
                    }
                },
                'genus': {
                    genus_name: {
                        'taxon_count': <int>,
                        'read_count': <int>,
                    }
                },
                'family': {
                    family_name: {
                        'taxon_count': <int>,
                        'read_count': <int>,
                    }
                }
            }
        }
    """
    # Initialize nested defaultdicts for kingdom groups and ranks
    kingdom_data = {
        group: {
            'species': defaultdict(
                lambda: {'taxon_count': 0, 'read_count': 0}
            ),
            'genus': defaultdict(
                lambda: {'taxon_count': 0, 'read_count': 0}
            ),
            'family': defaultdict(
                lambda: {'taxon_count': 0, 'read_count': 0}
            ),
        }
        for group in KINGDOM_GROUPS
    }

    with kraken_summary_path.open() as f:
        reader = csv.DictReader(f, delimiter='\t')

        for row in reader:
            full_lineage = row.get('full_lineage', '')
            full_lineage_ranks = row.get('full_lineage_ranks', '')
            reads = int(row.get('reads', 0))

            if not full_lineage or not full_lineage_ranks:
                continue

            # Split the lineage into taxonomic levels
            lineage_parts = [
                part.strip()
                for part in full_lineage.split(';')
                if part.strip()
            ]
            rank_parts = [
                part.strip().lower()
                for part in full_lineage_ranks.split(';')
                if part.strip()
            ]

            # Get taxonomic assignments
            rank_assignments = _assign_ranks_to_lineage(
                lineage_parts, rank_parts
            )

            # Determine kingdom group from lineage
            kingdom_group = _determine_kingdom_group(lineage_parts)

            if kingdom_group == 'na':
                kingdom_data[kingdom_group] = {
                    'read_count': reads,
                    'taxon_count': 1,
                }

            # Accumulate counts for species, genus, and family
            for rank in ['species', 'genus', 'family']:
                taxon = rank_assignments.get(rank)
                if taxon:
                    taxon_data = kingdom_data[kingdom_group][rank][taxon]
                    taxon_data['taxon_count'] += 1
                    taxon_data['read_count'] += reads

    # Convert to final structure with totals and top N for each rank
    result = {}

    for group in KINGDOM_GROUPS:
        group_result = {}

        # Process each rank (species, genus, family)
        for rank in ['species', 'genus', 'family']:
            sorted_taxa = dict(
                sorted(
                    kingdom_data[group][rank].items(),
                    key=lambda x: x[1]['read_count'],
                    reverse=True
                )
            )

            # Limit to TOP_N_TAXA and group the rest under 'other'
            top_taxa = dict(list(sorted_taxa.items())[:TOP_N_TAXA])
            other_taxa = dict(list(sorted_taxa.items())[TOP_N_TAXA:])
            other_reads = sum(
                taxon['read_count'] for taxon in other_taxa.values()
            )
            other_taxa_count = sum(
                taxon['taxon_count'] for taxon in other_taxa.values()
            )

            if other_reads > 0:
                top_taxa['other'] = {
                    'read_count': other_reads,
                    'taxon_count': other_taxa_count,
                }

            group_result[rank] = top_taxa

        # Calculate totals for the kingdom (based on species data)
        total_reads = sum(
            taxon['read_count']
            for taxon in group_result['species'].values()
        )
        total_taxa = sum(
            taxon['taxon_count']
            for taxon in group_result['species'].values()
        )

        result[group] = {
            'read_count': total_reads,
            'taxon_count': total_taxa,
            **group_result
        }

    na_read_percent = round(
        kingdom_data['na']['read_count']
        / sum([
            taxon['read_count']
            for taxon in result.values()
        ] + [kingdom_data['na']['read_count']])
        * 100, 2
    )
    result['na'] = {
        **kingdom_data['na'],
        'read_percent': na_read_percent,
    }

    return result


def _determine_kingdom_group(lineage_parts: list[str]) -> str:
    """Determine which kingdom group a lineage belongs to.

    Args:
        lineage_parts: List of taxonomic names from the lineage

    Returns:
        Kingdom group name: 'plant', 'animal', 'bacteria', 'fungi',
        'virus', 'other' or 'na' if not classified.
    """
    if lineage_parts[0] == 'NA':
        return 'na'

    higher_taxa = [
        part.lower() for part in lineage_parts
    ][:4]
    entity = higher_taxa[0]

    if entity == 'viruses':
        return 'virus'

    if len(higher_taxa) < 2:
        return 'other'

    domain = higher_taxa[1]

    if domain == 'bacteria':
        return 'bacteria'

    if domain == 'archaea':
        return 'other'

    # For eukaryotes, need to distinguish plant, animal, fungi
    if domain == 'eukaryota':
        if any(marker in higher_taxa for marker in PLANT_MARKERS):
            return 'plant'

        if any(marker in higher_taxa for marker in ANIMAL_MARKERS):
            return 'animal'

        if any(marker in higher_taxa for marker in FUNGI_MARKERS):
            return 'fungi'

    return 'other'


def _assign_ranks_to_lineage(
    lineage_parts: list[str],
    rank_parts: list[str],
) -> dict[str, str]:
    """Assign taxonomic ranks to lineage parts using rank labels.

    This function maps lineage parts to standard taxonomic ranks using
    the full_lineage_ranks data from the Kraken summary file.

    Args:
        lineage_parts: List of taxonomic names from the lineage
        rank_parts: List of rank labels corresponding to lineage_parts
                   (e.g., ['domain', 'kingdom', 'phylum', ..., 'species'])

    Returns:
        Dictionary mapping rank names to taxon names
    """
    assignments = {}

    if not lineage_parts or not rank_parts:
        return assignments

    if len(lineage_parts) != len(rank_parts):
        return assignments

    for taxon, rank in zip(lineage_parts, rank_parts):
        if rank in TAXONOMIC_RANKS:
            assignments[rank] = taxon

    return assignments


def get_top_taxa_by_rank(
    taxa_dict: dict[str, dict[str, dict[str, int]]],
    rank: str,
    top_n: int = 10,
) -> list[tuple[str, dict[str, int]]]:
    """Get the top N taxa for a given rank by read count.

    Args:
        taxa_dict: The nested dictionary returned by parse_kraken_taxonomy
        rank: The taxonomic rank to query (e.g., 'Species', 'Genus')
        top_n: Number of top taxa to return (default: 10)

    Returns:
        List of tuples (taxon_name, stats_dict) ordered by read_count

    Example:
        >>> taxa = parse_kraken_taxonomy('kraken_summary.txt')
        >>> top_species = get_top_taxa_by_rank(taxa, 'Species', top_n=5)
        >>> for name, stats in top_species:
        ...     print(f"{name}: {stats['read_count']} reads")
    """
    if rank not in taxa_dict:
        return []

    rank_data = taxa_dict[rank]
    sorted_items = sorted(
        rank_data.items(),
        key=lambda x: x[1]['read_count'],
        reverse=True
    )

    return sorted_items[:top_n]
