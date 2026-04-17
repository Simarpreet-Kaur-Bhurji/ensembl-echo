# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
generate_input_tsv.py
---------------------
Core input-processing library for ECHO pipeline step 1.

Responsibilities:
  1. load_species_info()         – parse the 4-column metadata TSV into a species map
  2. resolve_file_path()         – resolve relative/absolute FASTA paths
  3. combine_fastas_from_map()   – read per-species FASTAs, emit a single combined FASTA
  4. write_protein_metadata()    – join combined FASTA with species map → TSV + Parquet

Header format used throughout the pipeline:
    {protein_id}|{sps_name}_{taxon_id}_{gca}|{seq_len}

The composite species key (sps_name_taxon_id_gca) is the species identifier in:
  - the FASTA header (field 2, |-delimited)
  - the 'name' column of processed_input.parquet
  - the keys of the species_map dict returned by load_species_info()
"""

import os
import csv

import duckdb
import pandas as pd


# ---------------------------------------------------------------------------
# Low-level FASTA helpers
# ---------------------------------------------------------------------------


def write_sequence(out_f, header, seq, combined, line_width=60):
    """
    Write one sequence to the combined FASTA and append to the in-memory list.

    The sequence length is appended to the header so downstream code can read
    it back from the header without re-computing it:
        {protein_id}|{species_key}  →  {protein_id}|{species_key}|{seq_len}

    Args:
        out_f      : open file handle for the combined FASTA
        header     : header string (without leading '>')
        seq        : full sequence string
        combined   : list accumulating (header, seq) tuples in memory
        line_width : characters per sequence line (default 60; set via params.fasta_line_width)
    """
    seq_len = len(seq)
    header_with_len = f"{header}|{seq_len}"
    combined.append((header_with_len, seq))
    out_f.write(f">{header_with_len}\n")
    for i in range(0, seq_len, line_width):
        out_f.write(seq[i : i + line_width] + "\n")


def process_fasta_file(filepath, prefix, out_f, combined, line_width=60):
    """
    Stream-parse a single FASTA file (no index required), prepend the species
    prefix to each protein header, and write to the combined FASTA.

    Output header per sequence:
        {raw_protein_id}|{prefix}|{seq_len}

    Args:
        filepath   : path to the input FASTA
        prefix     : species key (sps_name_taxon_id_gca) prepended to headers
        out_f      : open file handle for the combined FASTA
        combined   : list accumulating (header, seq) tuples in memory
        line_width : characters per sequence line (passed through to write_sequence)
    """
    header = None
    seq_lines = []
    with open(filepath, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                # flush the previous record before starting a new one
                if header and seq_lines:
                    write_sequence(
                        out_f, header, "".join(seq_lines), combined, line_width
                    )
                # take only the first token so spaces in FASTA headers don't propagate
                raw_id = line[1:].strip().split()[0]
                header = f"{raw_id}|{prefix}"
                seq_lines = []
            else:
                seq_lines.append(line)
        # flush the final record in the file
        if header and seq_lines:
            write_sequence(out_f, header, "".join(seq_lines), combined, line_width)


# ---------------------------------------------------------------------------
# Metadata TSV loading
# ---------------------------------------------------------------------------


def load_species_info(tsv_file):
    """
    Load species metadata from the 4-column metadata TSV.

    Expected columns (tab-separated, with header row):
      - sps_name   : normalised species name (e.g. homo_sapiens)
      - taxon_id   : NCBI taxon ID
      - gca        : genome assembly accession, or NA
      - file_path  : path to the species FASTA — absolute, OR a bare filename /
                     relative path resolved against --input_fasta_dir at runtime

    Returns:
        dict keyed by the composite identifier  {sps_name}_{taxon_id}_{gca}.

    The composite key is used as the species prefix in every FASTA header and
    in the 'name' column of the output parquet, ensuring two assemblies of the
    same species (same sps_name, different gca) are always distinguishable.
    """
    mapping = {}
    print(f"[load_species_info] reading: {tsv_file}")
    with open(tsv_file, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        print(f"[load_species_info] columns detected: {reader.fieldnames}")
        for row in reader:
            sps_name = row["sps_name"].strip().replace(" ", "_")
            taxon_id = row["taxon_id"].strip()
            gca = row.get("gca", "NA").strip()
            # composite key uniquely identifies one assembly
            key = f"{sps_name}_{taxon_id}_{gca}"
            mapping[key] = {
                "name": key,  # stored as 'name' in the parquet
                "tax_id": taxon_id,
                "gca": gca,
                "file_path": row["file_path"].strip(),
            }
            print(
            f"  [load_species_info] loaded: key={key!r}  taxon_id={taxon_id}"
            f" file_path={mapping[key]['file_path']!r}"
            )
    print(f"[load_species_info] total species loaded: {len(mapping)}")
    return mapping


# ---------------------------------------------------------------------------
# File-path resolution
# ---------------------------------------------------------------------------


def resolve_file_path(file_path, input_fasta_dir=None):
    """
    Resolve a file_path entry from the metadata TSV to an accessible path.

    Resolution order:
      1. Absolute path in TSV  →  used as-is
      2. Relative path + input_fasta_dir given  →  joined with input_fasta_dir
      3. Relative path, no input_fasta_dir  →  relative to CWD (caller's responsibility)
    """
    if os.path.isabs(file_path):
        return file_path
    if input_fasta_dir:
        return os.path.join(input_fasta_dir, file_path)
    return file_path


# ---------------------------------------------------------------------------
# FASTA combining
# ---------------------------------------------------------------------------


def combine_fastas_from_map(
    species_map, combined_fasta_file, input_fasta_dir=None, line_width=60
):
    """
    Combine per-species FASTA files (listed in species_map) into one file.

    Each species' FASTA path comes from species_map[key]['file_path'], resolved
    via resolve_file_path().  The composite species key is embedded in every
    output header as the species prefix.

    Args:
        species_map         : dict returned by load_species_info()
        combined_fasta_file : output path for the merged FASTA
        input_fasta_dir     : optional base directory for relative file_path values
        line_width          : characters per sequence line in the output FASTA (default 60)

    Returns:
        list of (header, sequence) tuples (all sequences, in TSV order)
    """
    combined = []
    print(f"[combine_fastas_from_map] output: {combined_fasta_file}")
    print(
        f"[combine_fastas_from_map] input_fasta_dir (base for relative paths): {input_fasta_dir!r}"
    )
    print(f"[combine_fastas_from_map] fasta line width: {line_width}")

    with open(combined_fasta_file, "w", encoding="utf-8") as out_f:
        for species_key, info in species_map.items():
            filepath = resolve_file_path(info["file_path"], input_fasta_dir)
            print(
                f"  [combine_fastas_from_map] processing: {species_key!r}  ->  {filepath}"
            )
            if not os.path.exists(filepath):
                raise FileNotFoundError(
                    f"FASTA not found for {species_key!r}: {filepath}\n"
                    f"  (file_path in TSV: {info['file_path']!r}, input_fasta_dir: {input_fasta_dir!r})"
                )
            before = len(combined)
            process_fasta_file(filepath, species_key, out_f, combined, line_width)
            print(f"    -> {len(combined) - before} sequences added")

    print(f"[combine_fastas_from_map] total sequences: {len(combined)}")
    return combined


# ---------------------------------------------------------------------------
# Metadata output (TSV + Parquet)
# ---------------------------------------------------------------------------


def header_to_species_key(header):
    """
    Extract the composite species key from a combined FASTA header.

    Header format:  {protein_id}|{sps_name}_{taxon_id}_{gca}|{seq_len}
    Returns field [1], which is the key used to look up species_map.
    """
    try:
        return header.split("|")[1]
    except IndexError:
        return "unknown"


def write_protein_metadata(combined_fasta, species_map, output_tsv, _out_parquet):
    """
    Write per-protein metadata to a TSV and a Parquet file.

    For each sequence in combined_fasta:
      - protein_id and seq_len are parsed from the header
      - name, tax_id, gca are looked up from species_map via the composite key

    Output columns (both TSV and Parquet):
        protein_id, name, sequence_length, tax_id, gca

    The Parquet additionally stores the full header and sequence for downstream
    cluster-joining (parse_cluster_file uses header as the join key).

    Args:
        combined_fasta : list of (header, seq) tuples from combine_fastas_from_map()
        species_map    : dict from load_species_info()
        output_tsv     : path for the output TSV
        out_parquet    : path for the output Parquet (ignored here; derived from output_tsv)
    """
    rows = []
    unmatched = set()

    with open(output_tsv, "w", newline="", encoding="utf-8") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(
            [
                "protein_id",
                "name",
                "sequence_length",
                "tax_id",
                "gca",
            ]
        )

        for header, seq in combined_fasta:
            parts = header.split("|")
            if len(parts) < 3:
                raise ValueError(f"Header not in expected format: {header!r}")

            protein_id = parts[0]
            # seq_len was appended by write_sequence(); fall back to computing it
            seq_len_str = parts[-1]
            seq_len = int(seq_len_str) if seq_len_str.isdigit() else len(seq)

            # look up species metadata using the composite key embedded in the header
            species_key = header_to_species_key(header)
            info = species_map.get(species_key)
            if info is None:
                # warn once per unmatched key; write NA values so the run doesn't abort
                if species_key not in unmatched:
                    print(
                    f"  [write_protein_metadata] WARNING: no metadata match"
                    f" for key {species_key!r} — writing NA"
                    )
                    unmatched.add(species_key)
                info = {
                    "name": species_key,  # preserve the key so the header join still works
                    "tax_id": "NA",
                    "gca": "NA",
                }

            writer.writerow(
                [
                    protein_id,
                    info["name"],
                    seq_len,
                    info["tax_id"],
                    info["gca"],
                ]
            )

            # parquet row includes header + sequence for downstream cluster joining
            rows.append(
                [
                    header,
                    seq,
                    protein_id,
                    info["name"],
                    info["tax_id"],
                    info["gca"],
                    seq_len,
                ]
            )

    if unmatched:
        print(
        f"[write_protein_metadata] WARNING: {len(unmatched)}"
        f" unmatched species key(s): {sorted(unmatched)}"
        )

    df = pd.DataFrame(
        rows,
        columns=[
            "header",
            "sequence",
            "protein_id",
            "name",
            "tax_id",
            "gca",
            "seq_len",
        ],
    )

    # derive parquet path from TSV path (same stem, different extension)
    # use duckdb to write parquet to avoid pyarrow dependency
    parquet_file = output_tsv.rsplit(".", 1)[0] + ".parquet"
    con = duckdb.connect()
    con.register("df", df)
    con.execute(f"COPY df TO '{parquet_file}' (FORMAT PARQUET)")
    print(f"[write_protein_metadata] wrote {len(df)} proteins → {output_tsv}")
    print(f"[write_protein_metadata] wrote parquet               → {parquet_file}")
