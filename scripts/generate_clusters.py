#!/usr/bin/env python3

import argparse
import subprocess
import os
from collections import defaultdict
import glob
import json

def run_mmseqs(input_fasta, output_dir,
               min_seq_id=0.75, coverage=0.8, cov_mode=1, threads=16,
               singularity_image="/hps/nobackup/flicek/ensembl/compara/jitender/containers/mmseqs2_latest.sif"):
    """
    Run MMseqs2 easy-cluster using Singularity and return the cluster file path.
    """
    os.makedirs(output_dir, exist_ok=True)

    # MMseqs2 easy-cluster requires three arguments: input, output prefix, tmp
    output_prefix = os.path.join(output_dir, "mmseqs_results")
    tmp_dir = os.path.join(output_dir, "tmp")
    os.makedirs(tmp_dir, exist_ok=True)

    cmd = [
        "singularity", "exec",
        singularity_image,
        "mmseqs", "easy-cluster",
        input_fasta,
        output_prefix,
        tmp_dir,
        "--min-seq-id", str(min_seq_id),
        "-c", str(coverage),
        "--cov-mode", str(cov_mode),
        "--threads", str(threads)
    ]

    print("Running command:")
    print(" ".join(cmd))
    subprocess.run(cmd, check=True)
    print(f"MMseqs2 clustering completed. Results stored in: {output_dir}")

    # Return the full path to the cluster file
    cluster_file = f"{output_prefix}_cluster.tsv"
    return cluster_file


def parse_cluster_file(cluster_file):
    """
    Parse MMseqs2 easy-cluster tab-delimited output and return
    a dictionary in alfatclust style: keys '#Cluster N', values = list of protein IDs.

    Args:
        cluster_file (str): Path to the cluster file.

    Returns:
        dict: Dictionary where keys are '#Cluster N' and values are lists of protein IDs.
    """
    cluster_pool = defaultdict(list)

    # Read raw MMseqs2 cluster file
    with open(cluster_file, "r") as f:
        for line in f:
            cluster_id, protein = line.strip().split("\t")
            assert cluster_id, "No cluster ID found!"
            assert protein, "No protein ID found!"
            cluster_pool[cluster_id].append(protein)

    # Sanity check: ensure unique proteins in each cluster
    for cluster_id, members in cluster_pool.items():
        if len(members) != len(set(members)):
            raise ValueError(f"Duplicate protein IDs found in cluster {cluster_id}")

    # Convert to alfatclust style dictionary (#Cluster N)
    alfatclust_dict = {f"#Cluster {k}": members
                       for k, (_, members) in enumerate(cluster_pool.items(), start=1)}

    print(f"Parsed {len(alfatclust_dict)} clusters from {cluster_file}")

    json_output_path = os.path.join(os.path.dirname(cluster_file), "cluster_dict.json")
    with open(json_output_path, "w") as jf:
        json.dump(alfatclust_dict, jf, indent=2)
    print(f"Cluster dictionary saved as JSON: {json_output_path}")

    return alfatclust_dict
