import os
import pandas as pd
import json
import numpy as np
import plotly.express as px
import plotly.io as pio
import plotly.graph_objects as go


# ===============================
# Helper Functions
# ===============================


def load_cluster_parquet(cluster_parquet_path):
    """
    Load cluster data from a Parquet file and return a dictionary:
    {Cluster_ID: [header1, header2, ...]}
    """
    df = pd.read_parquet(cluster_parquet_path)
    clusters_dict = df.groupby("Cluster_ID")["header"].apply(list).to_dict()
    return clusters_dict


def load_tsv(tsv_path):
    """Load TSV into pandas DataFrame."""
    return pd.read_csv(tsv_path, sep="\t")


# ===============================
# Cluster Distribution Functions
# ===============================


def cluster_distribution_by_taxids(
    cluster_log_path, less_than_required_log_path, output_dir
):
    # Load data
    cluster_df = load_tsv(cluster_log_path)[["cluster_id", "num_unique_tax_ids"]]
    less_than_df = load_tsv(less_than_required_log_path)[["Cluster_ID", "Total_TaxIDs"]]

    # Rename column to match cluster_df
    less_than_df = less_than_df.rename(
        columns={"Total_TaxIDs": "num_unique_tax_ids", "Cluster_ID": "cluster_id"}
    )

    combined_df = pd.concat([cluster_df, less_than_df], ignore_index=True)
    combined_df["num_unique_tax_ids"] = combined_df["num_unique_tax_ids"].astype(int)

    # Group by num_unique_tax_ids
    distribution_df = (
        combined_df.groupby("num_unique_tax_ids")
        .size()
        .reset_index(name="Num_Clusters")
    )
    distribution_df = distribution_df.sort_values("num_unique_tax_ids")

    # Interactive line plot
    fig = px.line(
        distribution_df,
        x="num_unique_tax_ids",
        y="Num_Clusters",
        markers=True,
        labels={
            "num_unique_tax_ids": "Unique Tax IDs per cluster",
            "Num_Clusters": "Number of clusters",
        },
        title="Distribution of Clusters by Unique Tax IDs",
    )

    fig.update_traces(hovertemplate="Unique Tax IDs: %{x}<br>Clusters: %{y:.0f}")

    # Add label in top-right corner
    max_taxids = distribution_df["num_unique_tax_ids"].max()
    fig.add_annotation(
        text=f"Maximum number of unique tax IDs in a cluster = {max_taxids}",
        x=1,
        y=1,
        xref="paper",
        yref="paper",
        showarrow=False,
        font=dict(color="black", size=16),
        align="right",
        xanchor="right",
        yanchor="top",
    )

    # Save interactive HTML + static PNG
    output_html = os.path.join(output_dir, "cluster_distribution_by_taxids.html")
    output_png = os.path.join(output_dir, "cluster_distribution_by_taxids.png")
    fig.write_html(output_html)
    pio.write_image(fig, output_png, format="png", scale=2)

    print(f"Saved {output_html} and {output_png}")

    return distribution_df


def plot_proteins_per_cluster(cluster_parquet_path, output_dir):
    # Load clusters dict
    clusters_dict = load_cluster_parquet(cluster_parquet_path)

    # Compute number of proteins per cluster
    cluster_sizes = {k: len(v) for k, v in clusters_dict.items()}
    df = pd.DataFrame(
        list(cluster_sizes.items()), columns=["Cluster_ID", "Num_Proteins"]
    )

    # Exclude singletons
    df = df[df["Num_Proteins"] > 1]

    # Distribution
    dist = df.groupby("Num_Proteins").size().reset_index(name="Num_Clusters")
    dist = dist.sort_values("Num_Proteins")

    total_clusters = len(df)
    # Interactive line plot
    fig = px.line(
        dist,
        x="Num_Proteins",
        y=np.log10(dist["Num_Clusters"] + 1),
        markers=True,
        labels={
            "Num_Proteins": "Proteins per cluster",
            "y": "log10(Number of clusters + 1)",
        },
        title="Proteins per Cluster (Excluding Singletons)",
    )

    # Hover info with raw cluster counts
    fig.update_traces(
        hovertemplate="Proteins: %{x}<br>Clusters: %{text}", text=dist["Num_Clusters"]
    )

    # Add annotation (top-right corner)
    fig.add_annotation(
        text=f"Total clusters (excl. singletons) = {total_clusters}",
        x=1,
        y=1,
        xref="paper",
        yref="paper",
        showarrow=False,
        font=dict(color="black", size=14),
        align="right",
        xanchor="right",
        yanchor="top",
    )

    fig.update_layout(
        xaxis_title="Proteins per cluster",
        yaxis_title="Number of clusters",
        template="plotly_white",
    )
    # Save interactive HTML + static PNG
    output_html = os.path.join(output_dir, "cluster_distribution_by_proteins.html")
    output_png = os.path.join(output_dir, "cluster_distribution_by_proteins.png")
    fig.write_html(output_html)
    pio.write_image(fig, output_png, format="png", scale=2)

    print(f"Saved {output_html} and {output_png}")

    return dist


def cluster_summary(cluster_parquet_path, less_than_log, closest_rel_log, output_dir):
    clusters_dict = load_cluster_parquet(cluster_parquet_path)
    clusters_with_fewer_taxids = load_tsv(less_than_log)
    closest_relatives_df = load_tsv(closest_rel_log)

    total_clusters = len(clusters_dict.keys())
    total_singletons = sum(1 for v in clusters_dict.values() if len(v) == 1)
    num_fewer_taxid_clusters = len(clusters_with_fewer_taxids)
    num_closest_relatives_clusters = total_clusters - (
        total_singletons + num_fewer_taxid_clusters
    )

    # Write to file

    summary_file = os.path.join(output_dir, "cluster_summary.txt")
    mode = "a" if os.path.exists(summary_file) else "w"
    with open(summary_file, mode) as f:
        f.write("\n")
        f.write(
            f"Total number of clusters obatined from input fasta: {total_clusters}\n"
        )
        f.write(
            f"Singleton clusters (clusters with only one protein): {total_singletons}\n"
        )
        f.write(
            f"Number of clusters with fewer than the required tax IDs: {num_fewer_taxid_clusters}\n"
        )
        f.write(
            f"Closest relatives clusters (clusters with required tax IDs or more): {num_closest_relatives_clusters}\n"
        )

    print(f"\nCluster summary written to: {summary_file}")

    return {
        "total_clusters": total_clusters,
        "singleton_clusters": total_singletons,
        "clusters_with_fewer_taxids": num_fewer_taxid_clusters,
        "valid_clusters": num_closest_relatives_clusters,
    }


# ===============================
# ECHO Pipeline Stats Functions
# ===============================

import os


def get_total_sequences(fasta_file):
    """Count number of sequences in a FASTA file."""
    count = 0
    if os.path.exists(fasta_file):
        with open(fasta_file) as f:
            for line in f:
                if line.startswith(">"):
                    count += 1
    return count


def echo_pipeline_summary(output_dir, input_fasta_file):
    """
    Calculate ECHO pipeline sequence stats and write to echo_pipeline_summary.txt.

    Args:
        output_dir (str): Directory containing *_relatives.fa, discarded_singletons.fa, clusters_with_fewer_tax_ids.fa
        input_fasta_file (str): Original HCP input FASTA file.
    """
    total_input_seqs = get_total_sequences(input_fasta_file)

    # Find all *_relatives.fa files
    all_relatives_files = {
        fname.replace("_relatives.fa", ""): os.path.join(output_dir, fname)
        for fname in os.listdir(output_dir)
        if fname.endswith("_all_relatives.fa") and not fname.startswith(".")
    }

    # Count sequences per query
    retained_counts = {}
    for query, fasta_file in all_relatives_files.items():
        retained_counts[query] = get_total_sequences(fasta_file)

    # Count discarded singletons
    singletons_file = os.path.join(output_dir, "discarded_singletons.fa")
    singletons_count = get_total_sequences(singletons_file)

    # Count clusters with fewer tax IDs
    clusters_fewer_file = os.path.join(output_dir, ".clusters_with_fewer_tax_ids.fa")
    clusters_fewer_count = get_total_sequences(clusters_fewer_file)

    # All queris combined file number of sequences
    all_queries_combined_file = os.path.join(
        output_dir, "all_queries_relatives_combined.fa"
    )
    all_queries_combined_count = get_total_sequences(all_queries_combined_file)

    # Compute percentages
    retention_stats = {}
    for query, count in retained_counts.items():
        percent_retained = (count / total_input_seqs) * 100 if total_input_seqs else 0
        percent_discarded = 100 - percent_retained
        retention_stats[query] = {
            "Retained": count,
            "Percent Retained": percent_retained,
            "Percent Discarded": percent_discarded,
        }

    summary_stats = {
        "total_input": total_input_seqs,
        "retained_per_query": retained_counts,  # per-query sequences
        "singletons": singletons_count,  # common discarded
        "clusters_fewer_taxids": clusters_fewer_count,  # common fewer tax IDs
        "all_queries_combined": all_queries_combined_count,
    }

    # Write to echo_pipeline_summary.txt
    summary_file = os.path.join(output_dir, "echo_pipeline_summary.txt")
    with open(summary_file, "w") as f:
        f.write("ECHO Pipeline Summary\n")
        f.write("====================\n\n")
        f.write(f"Total sequences in input FASTA: {total_input_seqs}\n\n")
        f.write("Common files for all queries:\n")
        f.write(f"  Discarded singletons: {singletons_count} sequences\n")
        f.write(f"  Clusters with fewer tax IDs: {clusters_fewer_count} sequences\n\n")
        f.write("Sequences retained per query:\n")
        for query, stats in retention_stats.items():
            f.write(
                f"  {query} : {stats['Retained']} sequences "
                f"({stats['Percent Retained']:.2f}% retained, "
                f"{stats['Percent Discarded']:.2f}% discarded)\n"
            )
        f.write(
            "Note: Each query_all_relatives.fa file includes sequences from clusters_with_fewer_tax_ids.fa.\n\n"
        )
        f.write(
            f"Total sequences in all_queries_combined_relatives.fa: {all_queries_combined_count}\n"
        )
        f.write(
            "Note: This file combines all query_all_relatives.fa files and includes sequences from clusters_with_fewer_tax_ids.fa.\n"
        )

    print(f"ECHO pipeline summary written to {summary_file}")
    return summary_stats


def plot_relatives_sequence_counts(summary_stats, output_dir):
    """
    Plot sequences in *_all_relatives.fa per query along with discarded singletons.

    Args:
        summary_stats (dict): Output from echo_pipeline_summary function.
        output_dir (str): Directory to save plot.
    """
    # Extract data
    retained_counts = summary_stats["retained_per_query"]
    singletons_count = summary_stats["singletons"]

    # Categories and values
    categories = list(retained_counts.keys()) + ["Discarded singletons"]
    values = list(retained_counts.values()) + [singletons_count]

    # Colors: green for queries, red for singletons
    colors = ["green"] * len(retained_counts) + ["red"]

    # Create bar chart
    fig = go.Figure(
        go.Bar(
            x=categories,
            y=values,
            text=values,
            textposition="auto",
            marker_color=colors,
        )
    )

    fig.update_layout(
        title="Sequences per query and discarded singletons",
        xaxis_title="Queries / Category",
        yaxis_title="Number of sequences",
        template="plotly_white",
        height=500,
        width=950,
    )

    # Save
    png_path = os.path.join(output_dir, "relatives_per_query.png")
    fig.write_image(png_path)

    print(f"Plot saved as:{png_path}")


def get_diagnostic_stats_and_plots(input_fasta_file, output_dir):
    """
    Run all diagnostic stats and plots for ECHO pipeline.

    Args:
        input_fasta_file (str): Original HCP input FASTA.
        output_dir (str): Directory where outputs and intermediate files are stored.
    """

    # ------------------------------
    # Gather static files from output_dir
    # ------------------------------
    files = os.listdir(output_dir)

    cluster_parquet_path = (
        os.path.join(output_dir, "clusters.parquet")
        if "clusters.parquet" in files
        else None
    )
    closest_rel_log = (
        os.path.join(output_dir, "closest_relatives_log.tsv")
        if "closest_relatives_log.tsv" in files
        else None
    )
    less_than_log = (
        os.path.join(output_dir, "clusters_with_fewer_taxids_summary.tsv")
        if "clusters_with_fewer_taxids_summary.tsv" in files
        else None
    )
    singleton_log = (
        os.path.join(output_dir, "singleton_cluster_summary.tsv")
        if "singleton_cluster_summary.tsv" in files
        else None
    )

    # ------------------------------
    # Run diagnostics and plots
    # ------------------------------

    # Plot proteins per cluster (excluding singletons)
    plot_proteins_per_cluster(cluster_parquet_path, output_dir)

    # Distribution of clusters by unique tax IDs
    cluster_by_taxids = cluster_distribution_by_taxids(
        closest_rel_log, less_than_log, output_dir
    )

    # Cluster summary (print stats and write to file)ls

    cluster_summary(cluster_parquet_path, less_than_log, closest_rel_log, output_dir)

    # ECHO pipeline summary (calculate sequences retained, discarded, etc.)
    summary_stats = echo_pipeline_summary(output_dir, input_fasta_file)

    # Plot sequences in *_all_relatives.fa and discarded singletons
    plot_relatives_sequence_counts(summary_stats, output_dir)

    print("All diagnostic stats and plots generated successfully.")

    return {"cluster_by_taxids": cluster_by_taxids, "summary_stats": summary_stats}
