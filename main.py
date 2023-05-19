#!/usr/bin/env python3

import argparse
import logging
import sys
from datetime import datetime

from Path import PathObject

LOGGING_FORMAT = "%(levelname)-8s :: %(asctime)s :: %(message)s"
logging.basicConfig(
    format=LOGGING_FORMAT, filename=f"{datetime.now()}.log", level=logging.INFO
)


def get_args():
    parser = argparse.ArgumentParser(
        description="Path Object to process and store information from .fasta inputs for further analysis/plotting"
    )
    parser.add_argument(
        "-i",
        "--input",
        nargs="+",
        required=True,
        help="Input .fasta file(s)",
    )
    parser.add_argument(
        "-n",
        "--name",
        nargs="?",
        required=True,
        help="Name of path (for saving outputs, labeling plots)",
    )
    parser.add_argument(
        "-c",
        "--cdhit",
        nargs="?",
        required=True,
        help="CD-Hit .clstr output file path",
    )
    parser.add_argument(
        "-b",
        "--blast",
        nargs="?",
        required=True,
        help="Blastn output file using '-outfmt 7' formatting option",
    )

    return parser.parse_args()


def main():
    args = get_args()
    print(args.__dict__)
    logging.info(args)

    path_name = args.name
    path_files = args.input
    blast_file = args.blast
    cdhit_file = args.cdhit

    logging.info("path name: %s", path_name)
    logging.info("path files: %s", path_files)

    path_obj = PathObject(path_name, path_files, blast_file, cdhit_file)
    # info_df = path_obj.merge_tables()
    b_mask = path_obj.blast_df()
    # using class methods:
    # get kmer_df
    # get cd-hit-df
    # merge dfs directly in main

    print()


if __name__ == "__main__":
    sys.exit(main())
