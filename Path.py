import argparse
import itertools
import logging
from datetime import datetime
import subprocess
import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from sklearn.cluster import DBSCAN

LOGGING_FORMAT = "%(levelname)-8s :: %(asctime)s :: %(message)s"
logging.basicConfig(format=LOGGING_FORMAT, filename=f"{datetime.now()}.log", level=logging.INFO)


def get_args():
    parser = argparse.ArgumentParser(
        description="Path Object to process and store information from .fasta inputs for further analysis/plotting"
    )
    parser.add_argument(
        "-i",
        "--input",
        nargs='+',
        required=True,
        help="Input .fasta file(s)",
    )
    parser.add_argument(
        "-n",
        "--name",
        nargs='?',
        required=True,
        help="Name of path (for saving outputs, labeling plots)",
    )
    parser.add_argument(
        "-c",
        "--cd-hit",
        nargs='?',
        required=True,
        help="CD-Hit .clstr output file path"
    )
    parser.add_argument(
        "-b",
        "--blast",
        nargs='?',
        required=True,
        help="Blastn output file using '-outfmt 7' formatting option"
    )

    return parser.parse_args()


def call_cd_hit(input_file):
    logging.info("Calculating CD-Hit clusters...")
    outfile_path_out = subprocess.run(
        ["./cd-hit.sh", input_file],
        universal_newlines=True,
        capture_output=True,
        text=True,
    )
    # raise an error if it failed
    if outfile_path_out.stderr:
        logging.error(f'Error when calling cd-hit:')
        logging.error(outfile_path_out.stderr)
        raise RuntimeError(outfile_path_out.stderr)

    outfile_path = outfile_path_out.stdout

    return outfile_path


class PathObject:
    def __init__(self, path_name: str, input_files: list):
        logging.info('Initializing PathObject...')
        self.path_name = path_name
        self.input_files = input_files
        self.seq_dict = {
            "index": [],
            "input_file": [],
            "header": [],
            "sequence": [],
            "length": [],
        }

        for file in input_files:
            logging.info(f'Checking file: {file}')
            file_path = Path(file)
            with open(file, "r") as file_handle:
                for record in SeqIO.parse(file_handle, "fasta"):
                    self.seq_dict["input_file"].append(file_path.name)
                    self.seq_dict["index"].append(record.id)
                    self.seq_dict["header"].append(record.description)
                    self.seq_dict["sequence"].append(record.seq.__str__())
                    self.seq_dict["length"].append(len(record.seq))

        self.df = pd.DataFrame.from_dict(self.seq_dict)
        self.df.set_index("index", inplace=True)

        self.kmer_df = self._calculate_k_mer_frequencies()

        self.dbscan_df = self.cluster_dbscan()

        # self.cd_hit_files = []
        # for input_file in self.input_files:
        #     self.cd_hit_files.append(call_cd_hit(Path(input_file).resolve().__str__()))

    def _calculate_k_mer_frequencies(self, size=4):
        logging.info('calculating k-mer frequencies...')
        nucleotides = ["A", "T", "C", "G"]
        a = ["".join(x) for x in list(itertools.product(nucleotides, repeat=size))]
        kmer_df = pd.DataFrame(columns=a)

        for row in self.df.itertuples():
            logging.info(row.Index)
            kmer_dict = {x: 0 for x in a}

            for ii in range(0, row.length, 1):
                if ii + 4 > row.length:
                    break
                else:
                    s = row.sequence[ii : ii + size]
                    if "N" not in s:
                        kmer_dict[s] += 1

            for k, val in kmer_dict.items():
                kmer_dict[k] = [val / (row.length // size)]

            kmer_dict["index"] = row.Index
            new_row = pd.DataFrame.from_dict(kmer_dict)
            kmer_df = pd.concat([kmer_df, new_row])

        kmer_df.set_index("index", inplace=True)

        return kmer_df

    def cluster_dbscan(self):
        logging.info("Calculating DBSCAN clusters...")
        index = self.kmer_df.index
        dbscan_result = DBSCAN(eps=0.2, min_samples=1).fit(self.kmer_df.to_numpy())
        dbscan_clusters = dbscan_result.labels_
        dbscan_df = pd.DataFrame.from_dict(
            {"index": list(index), "dbscan_cluster": dbscan_clusters}, orient="columns"
        ).set_index("index")

        return dbscan_df

    def parse_cd_hit_clusters(self):
        logging.info("Parsing CD-Hit clusters...")
        print(self)
        pass

    def merge_tables(self):
        logging.info("Merging outputs...")
        kmer_df = self.kmer_df
        dbscan_df = self.dbscan_df
        info_df = self.df

        df = pd.concat([info_df, dbscan_df, kmer_df], axis=1)

        return df


def main():
    args = get_args()
    print(args.__dict__)
    logging.info(args)

    path_name = args.name
    path_files = args.input

    logging.info("path name: %s", path_name)
    logging.info("path files: %s", path_files)

    path_obj = PathObject(path_name, path_files)
    info_df = path_obj.merge_tables()
    # using class methods:
    # get kmer_df
    # get cd-hit-df
    # merge dfs directly in main

    print()


if __name__ == "__main__":
    sys.exit(main())
