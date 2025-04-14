"""
process_clusters_module.py

This module provides functionality to parse a cluster file and organize its data into a dictionary.

Functions:
- parse_cluster_file: Reads a cluster file and creates a dictionary where each cluster is represented 
  by a unique key, and its members are stored as a list.

Input File Format:
The cluster file should be a text file with the following structure:
1. Lines starting with "#Cluster <number>" indicate the start of a new cluster.
2. Subsequent lines contain the members of the cluster.

Example:
#Cluster 1
member1
member2
#Cluster 2
member3
member4
"""

def parse_cluster_file(cluster_file):
    """
    Parse the cluster file and create a dictionary of clusters.

    Args:
        cluster_file (str): Path to the cluster file.

    Returns:
        dict: A dictionary where the keys are cluster numbers and the values are lists of cluster members.
    """
    clusters = {}
    current_cluster = None

    with open(cluster_file, "r") as f:
        for line in f:
            line = line.strip()  # Remove leading/trailing whitespaces
            if line.startswith("#Cluster"):
                # Extract cluster number
                cluster_number = line.split()[1]
                current_cluster = cluster_number
                clusters.setdefault(cluster_number, set())
            elif current_cluster and line:
                # Append line to current cluster
                clusters[current_cluster].add(line)

    # Convert sets to lists before returning
    for cluster_number, values_set in clusters.items():
        clusters[cluster_number] = list(values_set)

    return clusters

