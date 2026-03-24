# analyses_config.py
import os

def file_exists(file_path):
    """Check if a file exists at the given path."""
    return os.path.isfile(file_path)

def is_file_empty(file_path):
    """Check if a file is empty."""
    return os.path.isfile(file_path) and os.path.getsize(file_path) == 0