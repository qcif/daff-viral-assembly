#!/usr/bin/env python

"""Build the HTML report from output data."""

import argparse

from report import report
from report.utils import existing_path


def main():
    """Parse the command line arguments and build the report."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--samplesheet",
        type=existing_path,
        help="Path to the samplesheet CSV file.",
    )
    
    parser.add_argument(
        "--default_params_file",
        type=existing_path,
        help=("Path to the default parameters YAML file."),
    )
    
    parser.add_argument(
        "--params_file",
        type=existing_path,
        help=("Path to the user parameters YAML file. If no params were"
              " modified, pass an empty file."),
    )
    parser.add_argument(
        '--versions',
        type=existing_path,
        help="The yml file containing the versions of all the tools used in the pipeline.",
    )
    parser.add_argument(
        "--analyst",
        help="The name of the analyst running the workflow.",
    )
    parser.add_argument(
        "--facility",
        help="The name of the facility where the workflow was submitted from.",
    )
    parser.add_argument(
        '--result_dir',
        type=existing_path,
        help="The directory containing the output data.",
    )

    args = parser.parse_args()
    report.render(
        args.result_dir,
        args.samplesheet,
        args.default_params_file,
        args.params_file,
        args.versions,
        args.analyst,
        args.facility,
    )


if __name__ == '__main__':
    main()
