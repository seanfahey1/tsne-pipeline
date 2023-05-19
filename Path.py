import itertools
import logging
import re
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from sklearn.cluster import DBSCAN


def merge_tables(df_list):
    logging.info("Merging outputs...")
    df = pd.concat([df_list], axis=1)

    return df


class PathObject:
    def __init__(self, path_name: str, input_files: list, blast_file: str, cdhit_cluster_file: str):
        logging.info("Initializing PathObject...")
        self.path_name = path_name
        self.input_files = input_files
        self.blast_file = blast_file
        self.cdhit_cluster_file = cdhit_cluster_file

        self.seq_dict = {
            "index": [],
            "input_file": [],
            "header": [],
            "sequence": [],
            "length": [],
        }

        for file in input_files:
            logging.info(f"Checking file: {file}")
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
        self.cd_hit_clusters = self.parse_cd_hit_clusters()
        self.blast_df = self.parse_blast_output()

        # self.cd_hit_files = []
        # for input_file in self.input_files:
        #     self.cd_hit_files.append(call_cd_hit(Path(input_file).resolve().__str__()))

    def _calculate_k_mer_frequencies(self, size=4):
        logging.info("calculating k-mer frequencies...")
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
                    s = row.sequence[ii: ii + size]
                    if "N" not in s:
                        kmer_dict[s] += 1

            for k, val in kmer_dict.items():
                kmer_dict[k] = [val / (row.length // size)]

            kmer_dict["index"] = row.Index
            new_row = pd.DataFrame.from_dict(kmer_dict)
            kmer_df = pd.concat([kmer_df, new_row])

        kmer_df.set_index("index", inplace=True)

        return kmer_df

    def dbscan_df(self):
        logging.info("Calculating DBSCAN clusters...")
        index = self.kmer_df.index
        dbscan_result = DBSCAN(eps=0.2, min_samples=1).fit(self.kmer_df.to_numpy())
        dbscan_clusters = dbscan_result.labels_
        dbscan_df = pd.DataFrame.from_dict(
            {"index": list(index), "dbscan_cluster": dbscan_clusters}, orient="columns"
        ).set_index("index", inplace=True)

        return dbscan_df

    def cdhit_df(self):
        logging.info("Parsing CD-Hit clusters...")
        with open(self.cdhit_cluster_file, 'r') as cluster_file:
            cluster_str = cluster_file.read()

        clusters = {}
        for i, cluster in enumerate(cluster_str.split(">Cluster ")[1:]):
            matches = re.findall(r'>\s*(.*?)\.\.\.', cluster)
            for match in matches:
                clusters[match] = i

        # TODO: fix this
        
        cdhit_df = pd.DataFrame.from_dict(clusters, orient="columns").set_index("contig", inplace=True)
        return clusters

    def blast_df(self):
        logging.info("Parsing blast results...")

        with open(self.blast_file, 'r') as blast_file:
            blast_str = blast_file.read()

        blast_results = {
            'contig': [],
            'hits': [],
            'phage': [],
            'human': [],
            'macrophage': [],
            'assignment': [],
            'blast_text': [],
        }
        sections = blast_str.split('\nQuery= ')
        for section in sections:
            section = section.strip()
            blast_results['contig'].append(section[:section.find(' ')])

            table_str = re.search(r'Value\s*(.*?)\s*\n>', section).group(1)
            table = [r.strip() for r in table_str.split('\n')]
            blast_results['blast_text'] = table_str

            phage = 0
            human = 0
            macrophage = 0

            for row in table:
                if 'macrophage' in row.lower():
                    macrophage += 1
                elif 'human' in row.lower() or 'homo sapien' in row.lower():
                    human += 1
                elif 'phage' in row.lower():
                    phage += 1
            blast_results['phage'].append(phage)
            blast_results['human'].append(human)
            blast_results['macrophage'].append(macrophage)
            blast_results['hits'].append(len(table))

            if human > 0 or macrophage > 0 :
                assignment = 'drop'
            elif phage > 0:
                assignment = 'phage'
            else:
                assignment = 'not-phage'
            blast_results['assignment'].append(assignment)

        blast_df = pd.DataFrame.from_dict(blast_results, orient="columns").set_index("contig", inplace=True)

        return blast_df

    def merge_tables(self):
        logging.info("Merging outputs...")
        kmer_df = self.kmer_df
        dbscan_df = self.dbscan_df
        info_df = self.df

        df = pd.concat([info_df, dbscan_df, kmer_df], axis=1)

        return df
