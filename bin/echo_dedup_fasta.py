#!/usr/bin/env python3
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
echo_dedup_fasta.py
-------------------
Merge one or more FASTA files and deduplicate by header.
Optionally also deduplicate by exact sequence content (--dedup_sequences).

Header dedup always runs — first occurrence kept, input order preserved.
Sequence dedup is off by default; when enabled, writes dedup_report.tsv
listing every dropped header, which header it duplicated, and sequence length.

The manifest is not affected — this is a FASTA-only post-processing step.
"""

import argparse


def dedup_fasta(in_paths, out_path, dedup_sequences=False):
    seen_headers = set()
    seen_seq_to_header = {}  # seq -> first header that claimed it
    total = header_dups = seq_dups = written = 0
    report_rows = []

    def flush(header, seq_lines):
        nonlocal total, header_dups, seq_dups, written
        total += 1
        seq = "".join(seq_lines)

        if header in seen_headers:
            header_dups += 1
            return

        if dedup_sequences and seq in seen_seq_to_header:
            seq_dups += 1
            seen_headers.add(header)
            report_rows.append((header, seen_seq_to_header[seq], len(seq)))
            return

        seen_headers.add(header)
        if dedup_sequences:
            seen_seq_to_header[seq] = header
        out.write(f">{header}\n" + "\n".join(seq_lines) + "\n")
        written += 1

    with open(out_path, "w", encoding="utf-8") as out:
        for in_path in in_paths:
            with open(in_path, encoding="utf-8") as fh:
                header = None
                seq_lines = []
                for line in fh:
                    line = line.rstrip("\n")
                    if line.startswith(">"):
                        if header is not None:
                            flush(header, seq_lines)
                        header = line[1:]
                        seq_lines = []
                    elif line:
                        seq_lines.append(line)
                if header is not None:
                    flush(header, seq_lines)

    print(f"[dedup_fasta] total: {total}  written: {written}  "
          f"header_dups_removed: {header_dups}  seq_dups_removed: {seq_dups}")

    if dedup_sequences and report_rows:
        with open("dedup_report.tsv", "w", encoding="utf-8") as rf:
            rf.write("dropped_header\tkept_header\tsequence_length\n")
            for dropped, kept, length in report_rows:
                rf.write(f"{dropped}\t{kept}\t{length}\n")
        print(f"[dedup_fasta] wrote dedup_report.tsv ({len(report_rows)} sequence duplicates)")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", nargs="+", required=True,
                        help="Input FASTA files, merged in order")
    parser.add_argument("--output", required=True,
                        help="Output deduplicated FASTA")
    parser.add_argument("--dedup_sequences", action="store_true",
                        help="Remove exact sequence duplicates (off by default)")
    args = parser.parse_args()
    dedup_fasta(args.inputs, args.output, dedup_sequences=args.dedup_sequences)


if __name__ == "__main__":
    main()
