# analyses_config.py
import csv
import os

VIRUS_VIROID_RENAME_MAP_PATH = os.path.join(
    os.path.dirname(__file__),
    "virus_viroid_name_map.tsv",
)


def file_exists(file_path):
    """Check if a file exists at the given path."""
    return os.path.isfile(file_path)

def is_file_empty(file_path):
    """Check if a file is empty."""
    return os.path.isfile(file_path) and os.path.getsize(file_path) == 0


def standardize_viruse_viroid_names(df, column):
    """Replace known virus and viroid names using the shared TSV map."""
    if column not in df.columns:
        return df

    with open(VIRUS_VIROID_RENAME_MAP_PATH, newline="", encoding="utf-8") as map_file:
        reader = csv.DictReader(map_file, delimiter="\t")
        species_column = (
            "species_name"
            if "species_name" in (reader.fieldnames or [])
            else "virus_species"
        )
        name_map = {
            row["virus_name"]: row[species_column]
            for row in reader
            if row.get("virus_name") and row.get(species_column)
        }

    standardized = df.copy()
    standardized[column] = standardized[column].replace(name_map)
    return standardized