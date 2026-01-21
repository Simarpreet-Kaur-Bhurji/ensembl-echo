#!/usr/bin/env python3
import argparse
import os
import math
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from datetime import datetime


# -----------------------------
# Helpers
# -----------------------------
def ensure_dir(path: str):
    os.makedirs(path, exist_ok=True)


def count_fasta_headers(path: str) -> int:
    n = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                n += 1
    return n


def read_tsv_rows(tsv_path: str) -> int:
    with open(tsv_path) as fh:
        next(fh, None)  # header
        return sum(1 for _ in fh)


def annotate_peak(ax, x_vals, y_vals, label_prefix="Peak"):
    """Annotate the max y point on a simple x/y series."""
    if len(y_vals) == 0:
        return
    i = int(np.argmax(y_vals))
    x = x_vals[i]
    y = y_vals[i]
    ax.scatter([x], [y], s=40)
    ax.annotate(
        f"{label_prefix}: x={x}, y={y}",
        xy=(x, y),
        xytext=(10, 10),
        textcoords="offset points",
        arrowprops=dict(arrowstyle="->", lw=1),
        fontsize=9,
    )


def savefig(fig, pdf_pages, png_path):
    fig.savefig(png_path, dpi=200)
    pdf_pages.savefig(fig)
    plt.close(fig)


def add_text_page(pdf_pages, title, content):
    """Add a text page to the PDF."""
    fig, ax = plt.subplots(figsize=(8.5, 11))  # Letter size
    ax.axis('off')
    ax.text(0.05, 0.95, title, fontsize=14, fontweight='bold', va='top')
    ax.text(0.05, 0.90, content, fontsize=10, va='top')
    pdf_pages.savefig(fig)
    plt.close(fig)


def html_img_block(title, png_name, one_liner, extra=""):
    return f"""
    <div class="card">
      <h2>{title}</h2>
      <p class="one">{one_liner}</p>
      {f"<p class='extra'>{extra}</p>" if extra else ""}
      <img src="{png_name}" alt="{title}">
    </div>
    """


# -----------------------------
# Plots
# -----------------------------
def plot_cluster_size_distribution(cluster_sizes: pd.Series, pdf_pages, png_path):
    """
    Plot distribution of proteins per cluster (excluding singletons).
    y axis is log10(count+1) for readability.
    """
    sizes = cluster_sizes[cluster_sizes > 1].astype(int)
    if sizes.empty:
        # still produce empty plot
        fig, ax = plt.subplots()
        ax.set_title("Cluster size distribution (proteins per cluster; excluding singletons)")
        ax.set_xlabel("Proteins per cluster")
        ax.set_ylabel("log10(#clusters + 1)")
        ax.text(0.5, 0.5, "No non-singleton clusters", ha="center", va="center", transform=ax.transAxes)
        savefig(fig, pdf_pages, png_path)
        return

    dist = sizes.value_counts().sort_index()
    x = dist.index.to_numpy()
    y = dist.values
    ylog = np.log10(y + 1)

    fig, ax = plt.subplots()
    ax.plot(x, ylog, marker="o", linestyle="-")
    ax.set_title("Cluster size distribution (proteins per cluster; excluding singletons)")
    ax.set_xlabel("Proteins per cluster")
    ax.set_ylabel("log10(#clusters + 1)")
    annotate_peak(ax, x, y, label_prefix="Peak (raw)")
    ax.grid(True, alpha=0.3)
    savefig(fig, pdf_pages, png_path)


def plot_unique_taxids_distribution(unique_taxids: pd.Series, pdf_pages, png_path):
    """
    Distribution of unique tax IDs per cluster.
    """
    ut = unique_taxids.astype(int)
    dist = ut.value_counts().sort_index()
    x = dist.index.to_numpy()
    y = dist.values

    fig, ax = plt.subplots()
    ax.plot(x, y, marker="o", linestyle="-")
    ax.set_title("Distribution of unique tax IDs per cluster")
    ax.set_xlabel("Unique tax IDs per cluster")
    ax.set_ylabel("#clusters")
    annotate_peak(ax, x, y, label_prefix="Peak")
    ax.grid(True, alpha=0.3)
    savefig(fig, pdf_pages, png_path)


def plot_retained_vs_discarded(total_clusters: int, singleton_clusters: int, few_taxid_clusters: int, valid_clusters: int, pdf_pages, png_path):
    """
    Simple 3-bar plot: singletons, fewer-taxid, valid clusters.
    """
    labels = ["Singleton clusters", "Fewer-taxid clusters", "Valid clusters"]
    values = [singleton_clusters, few_taxid_clusters, valid_clusters]

    fig, ax = plt.subplots()
    ax.bar(labels, values)
    ax.set_title("Cluster retention summary")
    ax.set_ylabel("#clusters")
    ax.set_xticklabels(labels, rotation=15, ha="right")

    # annotate counts above bars
    for i, v in enumerate(values):
        ax.text(i, v, str(v), ha="center", va="bottom", fontsize=9)

    ax.grid(True, axis="y", alpha=0.3)
    savefig(fig, pdf_pages, png_path)


def plot_taxonomic_distance_distribution(closest_log_tsv: str, pdf_pages, png_path):
    """
    Histogram of taxonomic distances from closest_relatives_log.tsv.
    """
    df = pd.read_csv(closest_log_tsv, sep="\t")
    if "closest_distances" not in df.columns and "distance" not in df.columns:
        # handle: either list-like column or pre-exploded distances
        fig, ax = plt.subplots()
        ax.set_title("Taxonomic distance distribution")
        ax.text(0.5, 0.5, "No distance column found in closest log", ha="center", va="center", transform=ax.transAxes)
        savefig(fig, pdf_pages, png_path)
        return

    # Your log stores lists in a column; try to parse that robustly
    distances = []
    if "distance" in df.columns:
        distances = df["distance"].dropna().tolist()
    else:
        # closest_distances column likely contains python-like lists as strings
        for v in df["closest_distances"].dropna().tolist():
            if isinstance(v, list):
                distances.extend(v)
            else:
                s = str(v).strip()
                # try to parse like "[1, 2, 3]"
                s = s.strip("[]")
                if s:
                    parts = [p.strip() for p in s.split(",")]
                    for p in parts:
                        try:
                            distances.append(float(p))
                        except ValueError:
                            pass

    distances = [d for d in distances if d is not None and not (isinstance(d, float) and math.isnan(d))]
    if not distances:
        fig, ax = plt.subplots()
        ax.set_title("Taxonomic distance distribution")
        ax.text(0.5, 0.5, "No distances found", ha="center", va="center", transform=ax.transAxes)
        savefig(fig, pdf_pages, png_path)
        return

    distances = np.array(distances, dtype=float)
    print(f"DEBUG: distances count={len(distances)}, min={np.min(distances)}, max={np.max(distances)}")
    fig, ax = plt.subplots()
    ax.hist(distances, bins=40)
    ax.set_title("Taxonomic distance distribution (closest relatives)")
    ax.set_xlabel("Taxonomic distance")
    ax.set_ylabel("#selected relatives")

    # peak bin annotation
    counts, bin_edges = np.histogram(distances, bins=40)
    peak_i = int(np.argmax(counts))
    peak_left = bin_edges[peak_i]
    peak_right = bin_edges[peak_i + 1]
    peak_count = int(counts[peak_i])
    ax.annotate(
        f"Peak bin: [{peak_left:.2f}, {peak_right:.2f}] → {peak_count}",
        xy=(0.7, 0.95),
        xycoords="axes fraction",
        fontsize=9,
        ha="left",
        va="top",
        bbox=dict(boxstyle="round,pad=0.3", alpha=0.15),
    )

    ax.grid(True, axis="y", alpha=0.3)
    savefig(fig, pdf_pages, png_path)


def plot_size_vs_taxid_scatter(cluster_sizes: pd.Series, unique_taxids: pd.Series, pdf_pages, png_path):
    """
    Scatter: proteins per cluster vs unique tax IDs.
    """
    x = cluster_sizes.astype(int).to_numpy()
    y = unique_taxids.astype(int).to_numpy()

    fig, ax = plt.subplots()
    ax.scatter(x, y, s=10, alpha=0.6)
    ax.set_title("Cluster size vs taxonomic diversity")
    ax.set_xlabel("Proteins per cluster")
    ax.set_ylabel("Unique tax IDs per cluster")
    ax.grid(True, alpha=0.3)

    # annotate a rough “peak density” using mode of binned sizes and taxids
    xb = np.clip(x, 1, 500)
    yb = np.clip(y, 1, 200)
    hx, xedges, yedges = np.histogram2d(xb, yb, bins=[50, 40])
    pi = np.unravel_index(np.argmax(hx), hx.shape)
    x0 = (xedges[pi[0]] + xedges[pi[0]+1]) / 2
    y0 = (yedges[pi[1]] + yedges[pi[1]+1]) / 2
    ax.scatter([x0], [y0], s=60)
    ax.annotate(
        f"Highest density ~ ({x0:.0f}, {y0:.0f})",
        xy=(x0, y0),
        xytext=(10, 10),
        textcoords="offset points",
        arrowprops=dict(arrowstyle="->", lw=1),
        fontsize=9,
    )

    savefig(fig, pdf_pages, png_path)


# -----------------------------
# HTML report
# -----------------------------
def write_html(out_html: str, title: str, run_meta: dict, cards_html: str):
    meta_rows = "\n".join([f"<tr><th>{k}</th><td>{v}</td></tr>" for k, v in run_meta.items()])
    html = f"""<!doctype html>
<html>
<head>
<meta charset="utf-8">
<title>{title}</title>
<style>
  body {{ font-family: Arial, sans-serif; margin: 24px; }}
  h1 {{ margin-bottom: 6px; }}
  .muted {{ color: #555; }}
  table {{ border-collapse: collapse; margin-top: 10px; }}
  th, td {{ border: 1px solid #ddd; padding: 8px; }}
  th {{ text-align: left; background: #f6f6f6; }}
  .grid {{ display: grid; grid-template-columns: 1fr; gap: 18px; margin-top: 18px; }}
  .card {{ border: 1px solid #ddd; border-radius: 10px; padding: 14px; }}
  .one {{ font-weight: 600; margin-top: 6px; }}
  .extra {{ margin-top: 8px; }}
  img {{ width: 100%; max-width: 1100px; height: auto; border: 1px solid #eee; border-radius: 6px; margin-top: 12px; }}
</style>
</head>
<body>
  <h1>{title}</h1>
  <p class="muted">Generated: {datetime.now().isoformat(timespec="seconds")}</p>

  <h2>Run summary</h2>
  <table>
    {meta_rows}
  </table>

  <div class="grid">
    {cards_html}
  </div>
</body>
</html>
"""
    with open(out_html, "w") as f:
        f.write(html)


# -----------------------------
# Main
# -----------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--clusters_parquet", required=True)
    ap.add_argument("--remaining_clusters_parquet", required=True)
    ap.add_argument("--singleton_summary_tsv", required=True)
    ap.add_argument("--few_taxids_summary_tsv", required=True)
    ap.add_argument("--closest_log_tsv", required=True)
    ap.add_argument("--input_fasta", required=True)
    ap.add_argument("--cluster_summary_txt")
    ap.add_argument("--pipeline_summary_txt")
    ap.add_argument("--with_singletons", action="store_true")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    ensure_dir(args.outdir)

    # Load clusters
    cdf = pd.read_parquet(args.clusters_parquet)
    total_clusters = int(cdf["Cluster_ID"].nunique())
    cluster_sizes = cdf.groupby("Cluster_ID").size()

    # unique tax ids per cluster
    unique_taxids = cdf.groupby("Cluster_ID")["tax_id"].nunique()

    # Load remaining clusters to get actual valid clusters
    rdf = pd.read_parquet(args.remaining_clusters_parquet)
    valid_clusters = int(rdf["Cluster_ID"].nunique())

    # cluster discard counts from summary TSVs (cluster counts)
    singleton_clusters = read_tsv_rows(args.singleton_summary_tsv)  # number of singleton clusters
    few_taxid_clusters = read_tsv_rows(args.few_taxids_summary_tsv)  # number of few-taxid clusters

    # seq counts
    total_input_seqs = count_fasta_headers(args.input_fasta)

    # Create PDF for plots
    pdf_path = os.path.join(args.outdir, "diagnostics.pdf")
    with PdfPages(pdf_path) as pdf_pages:
        # Add summary page if files provided
        summary_content = ""
        if args.cluster_summary_txt and os.path.exists(args.cluster_summary_txt):
            with open(args.cluster_summary_txt, 'r') as f:
                summary_content += f.read() + "\n\n"
        if args.pipeline_summary_txt and os.path.exists(args.pipeline_summary_txt):
            with open(args.pipeline_summary_txt, 'r') as f:
                summary_content += f.read() + "\n\n"
        
        if summary_content:
            plot_descriptions = """
Plot Descriptions:

1. Cluster Size Distribution: Shows the distribution of proteins per cluster, excluding singletons. 
   The y-axis is log-transformed for better readability. Peaks indicate common cluster sizes.

2. Unique Tax IDs Distribution: Displays how many unique taxonomic IDs are present in each cluster.
   Helps understand taxonomic diversity within clusters.

3. Cluster Retention Summary: Bar chart showing total clusters, discarded clusters (singletons and few-taxid), 
   and retained valid clusters.

4. Taxonomic Distance Distribution: Histogram of taxonomic distances for closest relatives.
   Lower distances indicate closer taxonomic relationships.

5. Size vs Tax ID Scatter: Scatter plot correlating cluster size with number of unique taxonomic IDs.
   Can reveal patterns between cluster size and taxonomic diversity.

6. (If applicable) Additional plots for singletons or other metrics.
"""
            full_content = summary_content + plot_descriptions
            add_text_page(pdf_pages, "ECHO Pipeline Diagnostics Summary", full_content)
        
        # Make plots
        plot_cluster_size_distribution(cluster_sizes, pdf_pages, os.path.join(args.outdir, "plot1_cluster_size_distribution.png"))
        plot_unique_taxids_distribution(unique_taxids, pdf_pages, os.path.join(args.outdir, "plot2_unique_taxids_distribution.png"))
        plot_retained_vs_discarded(total_clusters, singleton_clusters, few_taxid_clusters, valid_clusters, pdf_pages, os.path.join(args.outdir, "plot3_cluster_retention_summary.png"))
        plot_taxonomic_distance_distribution(args.closest_log_tsv, pdf_pages, os.path.join(args.outdir, "plot5_taxonomic_distance_distribution.png"))
        plot_size_vs_taxid_scatter(cluster_sizes, unique_taxids, pdf_pages, os.path.join(args.outdir, "plot6_size_vs_taxid_scatter.png"))

    # Run meta table
    run_meta = {
        "Total input sequences (FASTA headers)": total_input_seqs,
        "Total clusters": total_clusters,
        "Singleton clusters": singleton_clusters,
        "Fewer-taxid clusters": few_taxid_clusters,
        "Valid clusters": valid_clusters,
        "with_singletons": args.with_singletons,
        "clusters.parquet": os.path.basename(args.clusters_parquet),
        "remaining_clusters.parquet": os.path.basename(args.remaining_clusters_parquet),
        "closest_relatives_log.tsv": os.path.basename(args.closest_log_tsv),
    }

    # Also write a tiny text summary for convenience
    with open(os.path.join(args.outdir, "diagnostics_summary.txt"), "w") as f:
        for k, v in run_meta.items():
            f.write(f"{k}\t{v}\n")

    print(f"[OK] Wrote diagnostics PDF and summary to: {args.outdir}")


if __name__ == "__main__":
    main()

