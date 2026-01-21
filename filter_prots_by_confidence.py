#!/usr/bin/env python3
import argparse
from pathlib import Path
import logging
import sys

def setup_logger(log_path):
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(log_path),
            logging.StreamHandler(sys.stdout),
        ],
    )


def load_allowed_gene_ids(conf_csv: Path) -> set[str]:
    allowed = set()
    with conf_csv.open() as fh:
        for line in fh:
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 4:
                continue
            gene_id = cols[3]              # col 3 (0-based) = gene_id
            confidence = cols[-1].lower()  # last column like "ensembl low"
            if confidence in {"high", "medium"}:
                allowed.add(gene_id)
    return allowed


def parse_gff_attr(attr: str) -> dict:
    d = {}
    for part in attr.strip().split(";"):
        if not part:
            continue
        if "=" in part:
            k, v = part.split("=", 1)
            d[k] = v
    return d


def load_transcript_to_gene_map(gff_path: Path) -> dict[str, str]:
    """
    From mRNA lines:
      ... mRNA ... Parent=<gene_id>;ID=<transcript_id>
    return {transcript_id: gene_id}
    """
    t2g = {}
    with gff_path.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            feature = cols[2]
            if feature != "mRNA":
                continue
            attrs = parse_gff_attr(cols[8])
            tid = attrs.get("ID")
            gid = attrs.get("Parent")
            if tid and gid:
                t2g[tid] = gid
    return t2g


def transcript_id_from_fasta_header(header_line: str) -> str:
    """
    Header line includes leading '>'.
    Your fasta looks like:
      >ENSZLMT00000000006|...
    so transcript is first token before '|'.
    """
    h = header_line[1:].strip()
    return h.split("|", 1)[0]


def filter_fasta_by_gene_conf(
    prots_fa: Path,
    allowed_genes: set[str],
    t2g: dict[str, str],
    out_fa: Path,
) -> tuple[int, int, int]:
    total = 0
    kept = 0
    unmapped = 0

    with prots_fa.open() as inp, out_fa.open("w") as out:
        write_seq = False
        for line in inp:
            if line.startswith(">"):
                total += 1
                tid = transcript_id_from_fasta_header(line)
                gid = t2g.get(tid)
                if gid is None:
                    unmapped += 1
                    write_seq = False
                    continue
                if gid in allowed_genes:
                    kept += 1
                    write_seq = True
                    out.write(line)
                else:
                    write_seq = False
            else:
                if write_seq:
                    out.write(line)

    return total, kept, unmapped


def find_by_prefix(files: list[Path], prefix: str, suffix: str) -> Path | None:
    candidates = [p for p in files if p.name.startswith(prefix) and p.name.endswith(suffix)]
    if not candidates:
        return None
    return sorted(candidates)[0]


def main():
    ap = argparse.ArgumentParser(
        description="Filter *.prots.fa by gene confidence using confidence.csv + transcript->gene mapping from GFF."
    )
    ap.add_argument("--input_dir", required=True, help="Directory containing *.prots.fa, *_confidence.csv, *_genome_annotations.gff")
    ap.add_argument("--out_dir", help="Output directory (default: input_dir)")
    args = ap.parse_args()

    in_dir = Path(args.input_dir)
    out_dir = Path(args.out_dir) if args.out_dir else in_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    prots_files = sorted(in_dir.glob("*.prots.fa"))
    conf_files = sorted(in_dir.glob("*_confidence.csv"))
    gff_files = sorted(in_dir.glob("*_genome_annotations.gff"))

    if not prots_files:
        raise SystemExit(f"[ERROR] No *.prots.fa found in {in_dir}")

    n_ok = 0
    n_skip = 0

    log_file = out_dir / "filter_prots_by_confidence.log"
    setup_logger(log_file)

    for prots_fa in prots_files:
        prefix = prots_fa.name.replace(".prots.fa", "")

        conf_csv = find_by_prefix(conf_files, prefix, "_confidence.csv")
        gff = find_by_prefix(gff_files, prefix, "_genome_annotations.gff")

        if conf_csv is None:
            print(f"[WARN] Missing confidence.csv for {prots_fa.name} (prefix={prefix}) — skipping")
            n_skip += 1
            continue
        if gff is None:
            print(f"[WARN] Missing genome_annotations.gff for {prots_fa.name} (prefix={prefix}) — skipping")
            n_skip += 1
            continue

        allowed_genes = load_allowed_gene_ids(conf_csv)
        t2g = load_transcript_to_gene_map(gff)

        out_fa = out_dir / prots_fa.name.replace(".prots.fa", ".prots.hm.fa")
        total, kept, unmapped = filter_fasta_by_gene_conf(prots_fa, allowed_genes, t2g, out_fa)


        logging.info(
            f"[OK] {prots_fa.name} -> {out_fa.name} | kept {kept}/{total} "
            f"(unmapped transcripts: {unmapped}, allowed genes: {len(allowed_genes)}, t2g size: {len(t2g)})"
        )
        n_ok += 1

    print(f"[DONE] Completed: {n_ok} written, {n_skip} skipped")


if __name__ == "__main__":
    main()
