"""Parse taxonomic classifications from Kraken summary files."""

import csv
from collections import defaultdict
from pathlib import Path
from typing import Union


# Define the standard taxonomic ranks in order
TAXONOMIC_RANKS = [
    'Domain',
    'Kingdom',
    'Phylum',
    'Class',
    'Order',
    'Family',
    'Genus',
    'Species',
]


def parse_kraken_taxonomy(
    kraken_summary_path: Union[str, Path],
) -> dict[str, dict[str, dict[str, int]]]:
    """Parse Kraken summary file and enumerate taxa by rank.

    Args:
        kraken_summary_path: Path to the Kraken summary TSV file

    Returns:
        A nested dictionary structure:
        {
            rank: {
                taxon: {
                    'taxon_count': <int>,  # number of rows with this taxon
                    'read_count': <int>,   # sum of reads for this taxon
                }
            }
        }
        Taxa within each rank are ordered by descending read abundance.

    Example:
        >>> result = parse_kraken_taxonomy('kraken_summary.txt')
        >>> result['Domain']['Eukaryota']['read_count']
        12345678
    """
    kraken_summary_path = Path(kraken_summary_path)

    # Initialize nested defaultdicts to accumulate data
    taxa_data = {rank: defaultdict(lambda: {'taxon_count': 0, 'read_count': 0})
                 for rank in TAXONOMIC_RANKS}

    with kraken_summary_path.open() as f:
        reader = csv.DictReader(f, delimiter='\t')

        for row in reader:
            full_lineage = row.get('full_lineage', '')
            reads = int(row.get('reads', 0))

            if not full_lineage:
                continue

            # Split the lineage into taxonomic levels
            lineage_parts = [
                part.strip()
                for part in full_lineage.split(';')
                if part.strip()
            ]

            # Map lineage parts to ranks
            # The first part is typically "cellular organisms" or root
            # followed by Domain, then standard ranks
            rank_assignments = _assign_ranks_to_lineage(lineage_parts)

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


def _assign_ranks_to_lineage(lineage_parts: list[str]) -> dict[str, str]:
    """Assign taxonomic ranks to lineage parts.

    This function attempts to map lineage parts to standard ranks.
    The mapping is based on position and known patterns in Kraken output.

    Common Kraken lineage pattern:
    0: cellular organisms (root - skip)
    1: Domain (Eukaryota, Bacteria, Archaea, Viruses)
    2+: Variable intermediate ranks, then standard ranks

    Args:
        lineage_parts: List of taxonomic names from the lineage

    Returns:
        Dictionary mapping rank names to taxon names
    """
    assignments = {}

    if not lineage_parts:
        return assignments

    # Skip "cellular organisms" if it's the first element
    start_idx = 1 if lineage_parts[0] == 'cellular organisms' else 0

    if len(lineage_parts) <= start_idx:
        return assignments

    # Position 1 (or 0 if no "cellular organisms") is typically Domain
    domain_idx = start_idx
    assignments['Domain'] = lineage_parts[domain_idx]

    # Try to identify standard ranks from the end of the lineage
    # Species is typically the last element
    if len(lineage_parts) > domain_idx + 1:
        assignments['Species'] = lineage_parts[-1]

    # Genus is typically second-to-last (if we have at least 2 elements)
    if len(lineage_parts) > domain_idx + 2:
        assignments['Genus'] = lineage_parts[-2]

    # For intermediate ranks, we use heuristics based on position
    # This is approximate since Kraken lineages can vary
    remaining_parts = lineage_parts[domain_idx + 1:-2]
    intermediate_ranks = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family']

    # Map remaining parts to intermediate ranks
    # taking the first N parts for the first N intermediate ranks
    for i, part in enumerate(remaining_parts):
        if i < len(intermediate_ranks):
            assignments[intermediate_ranks[i]] = part

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
