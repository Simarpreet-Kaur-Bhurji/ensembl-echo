"""
test_dedup.py
-------------
Unit tests for echo_dedup_fasta.py.

Tests:
  - Header deduplication always runs
  - Sequence deduplication only when --dedup_sequences is passed
  - dedup_report.tsv written only when --dedup_sequences is active and dups exist
  - Multi-file merge deduplicates across files

Headers use the real pipeline format: {prefix}|{species}_{tax_id}_{assembly}|{seq_len}

Run with:
  pytest test/test_dedup.py -v
"""

import logging
import os
import subprocess
import sys

import pytest

log = logging.getLogger("echo.test.dedup")

SCRIPT = os.path.join(os.path.dirname(__file__), "..", "bin", "echo_dedup_fasta.py")

# Real-format headers from the pipeline (format: prefix|species_taxid_assembly|length)
HDR_YAR = "yar02|yarrowia_lipolytica_4952_GCA_000002525.1|60"
HDR_CAN = "can02|candida_albicans_5476_GCA_000182965.3|56"
HDR_SCH = "sch02|schizosaccharomyces_pombe_4896_GCA_000002945.2|61"
HDR_ASP = "asp02|aspergillus_niger_5061_GCA_000002655.1|61"
HDR_NEU = "neu02|neurospora_crassa_5141_GCA_000182925.2|55"
HDR_SING = "asp05|aspergillus_niger_5061_GCA_000002655.1|52"

SEQ_A = "ACDEFGHIKLMNPQRSTVWY"
SEQ_B = "MKLVPQRSTEWYNACDFGHI"
SEQ_C = "GHIKLMNPQRSTVWYACDFE"

# Real sequences from the pipeline (from nextflow_test/saccharomyces_cerevisiae_all_relatives.fa)
# Used to test dedup with realistic protein data
REAL_SEQ_YAR02 = "QNQCWSWEEDDIAMSGSHRTHIYGRVDGCWSDPCSCQCWHTAACFTRGSLKPQLHKEEGA"
REAL_SEQ_CAN02 = "AEQCWSWEEDDIAMHGSWRSHMYGRVDGCTSDPKSCQCWHTAACFTRGPRKPQLHL"
REAL_SEQ_SCH02 = "QEQCWTWEEDDIAMSGSWRTHSYGRVDGCTSSPCSCCCWYTAACFTRGPRKPQLHLLKGDG"
REAL_SEQ_ASP02 = "MEQCWSWREDDIAMSGSWRTHIYDRVDGCTSDPCSCQDWHTAACVTRGPRKPQLHVPEKVS"
REAL_SEQ_NEU02 = "FEQCWSWEEDIIAMSGSWRTHIYGRMDGETSDPCSCQCWHTAACFTRGPRKLQLH"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def write_fasta(path, records):
    """records: list of (header, sequence)"""
    with open(path, "w", encoding="utf-8") as fh:
        for header, seq in records:
            fh.write(f">{header}\n{seq}\n")


def read_fasta_headers(path):
    with open(path, encoding="utf-8") as fh:
        return [line[1:].strip() for line in fh if line.startswith(">")]


def read_fasta_sequences(path):
    headers, seqs, cur_seq = [], [], []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if cur_seq:
                    seqs.append("".join(cur_seq))
                    cur_seq = []
                headers.append(line[1:])
            elif line:
                cur_seq.append(line)
    if cur_seq:
        seqs.append("".join(cur_seq))
    return dict(zip(headers, seqs))


def run_dedup(tmp_path, inputs, extra_args=None):
    out = tmp_path / "out.fa"
    cmd = [sys.executable, SCRIPT, "--inputs"] + [str(i) for i in inputs] + ["--output", str(out)]
    if extra_args:
        cmd += extra_args
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(tmp_path), check=False)
    assert result.returncode == 0, f"Script failed:\n{result.stderr}\n{result.stdout}"
    return out, result.stdout


# ---------------------------------------------------------------------------
# Header dedup
# ---------------------------------------------------------------------------

def test_duplicate_headers_removed(tmp_path):
    """Same header appearing twice → only first occurrence written."""
    fa = tmp_path / "input.fa"
    write_fasta(fa, [
        (HDR_YAR, SEQ_A),
        (HDR_CAN, SEQ_B),
        (HDR_YAR, SEQ_C),   # duplicate header, different sequence
    ])
    log.debug("input: 3 records, HDR_YAR appears twice with different sequences")
    out, stdout = run_dedup(tmp_path, [fa])
    headers = read_fasta_headers(out)
    log.debug("output headers (%d): %s", len(headers), headers)
    log.debug("dedup stdout: %s", stdout.strip())
    assert headers.count(HDR_YAR) == 1
    assert len(headers) == 2
    # First occurrence kept
    seqs = read_fasta_sequences(out)
    assert seqs[HDR_YAR] == SEQ_A
    assert "header_dups_removed: 1" in stdout
    log.info("PASS: duplicate header removed, first occurrence kept (%s)", HDR_YAR)


def test_unique_headers_all_retained(tmp_path):
    """No duplicate headers → all sequences written unchanged."""
    fa = tmp_path / "input.fa"
    records = [
        (HDR_YAR, SEQ_A),
        (HDR_CAN, SEQ_B),
        (HDR_SCH, SEQ_C),
        (HDR_ASP, SEQ_A),
        (HDR_NEU, SEQ_B),
    ]
    write_fasta(fa, records)
    log.debug("input: %d unique-header records", len(records))
    out, stdout = run_dedup(tmp_path, [fa])
    headers = read_fasta_headers(out)
    log.debug("output headers (%d): %s", len(headers), headers)
    assert headers == [h for h, _ in records]
    assert "header_dups_removed: 0" in stdout
    log.info("PASS: all %d unique headers retained", len(records))


# ---------------------------------------------------------------------------
# Sequence dedup (opt-in)
# ---------------------------------------------------------------------------

def test_sequence_dedup_off_by_default(tmp_path):
    """Same sequence, different headers → both retained when --dedup_sequences not set."""
    fa = tmp_path / "input.fa"
    write_fasta(fa, [
        (HDR_YAR, SEQ_A),
        (HDR_CAN, SEQ_A),   # same sequence, different header
    ])
    log.debug("input: 2 records sharing identical sequence SEQ_A; no --dedup_sequences flag")
    out, stdout = run_dedup(tmp_path, [fa])
    headers = read_fasta_headers(out)
    log.debug("output headers (%d): %s | stdout: %s", len(headers), headers, stdout.strip())
    assert len(headers) == 2
    assert "seq_dups_removed: 0" in stdout
    log.info("PASS: seq dedup off by default — both records retained")


def test_sequence_dedup_removes_exact_duplicates(tmp_path):
    """Same sequence, different headers → second dropped when --dedup_sequences."""
    fa = tmp_path / "input.fa"
    write_fasta(fa, [
        (HDR_YAR, SEQ_A),
        (HDR_CAN, SEQ_A),   # exact duplicate of HDR_YAR
        (HDR_SCH, SEQ_B),   # unique
    ])
    log.debug("input: YAR+CAN share SEQ_A, SCH is unique — running with --dedup_sequences")
    out, stdout = run_dedup(tmp_path, [fa], extra_args=["--dedup_sequences"])
    headers = read_fasta_headers(out)
    log.debug("output headers (%d): %s | stdout: %s", len(headers), headers, stdout.strip())
    assert len(headers) == 2
    assert HDR_YAR in headers
    assert HDR_CAN not in headers   # dropped — HDR_YAR was first
    assert HDR_SCH in headers
    assert "seq_dups_removed: 1" in stdout
    log.info("PASS: seq dup removed — kept=%s dropped=%s", HDR_YAR, HDR_CAN)


def test_sequence_dedup_writes_report(tmp_path):
    """dedup_report.tsv is written when --dedup_sequences is active and dups exist."""
    import pandas as pd
    fa = tmp_path / "input.fa"
    write_fasta(fa, [
        (HDR_YAR, SEQ_A),
        (HDR_CAN, SEQ_A),   # duplicate
    ])
    log.debug("input: YAR+CAN share SEQ_A (len=%d); expecting dedup_report.tsv", len(SEQ_A))
    run_dedup(tmp_path, [fa], extra_args=["--dedup_sequences"])
    report = tmp_path / "dedup_report.tsv"
    assert report.exists(), "dedup_report.tsv should be created when seq dups exist"
    df = pd.read_csv(report, sep="\t")
    log.debug("report columns: %s | rows: %d", list(df.columns), len(df))
    log.debug("report row: dropped=%s kept=%s len=%s",
              df.iloc[0]["dropped_header"], df.iloc[0]["kept_header"], df.iloc[0]["sequence_length"])
    assert list(df.columns) == ["dropped_header", "kept_header", "sequence_length"]
    assert len(df) == 1
    assert df.iloc[0]["dropped_header"] == HDR_CAN
    assert df.iloc[0]["kept_header"]    == HDR_YAR
    assert df.iloc[0]["sequence_length"] == len(SEQ_A)
    log.info("PASS: dedup_report.tsv written with correct dropped/kept/length")


def test_no_report_when_no_seq_dups(tmp_path):
    """dedup_report.tsv is NOT written when there are no sequence duplicates."""
    fa = tmp_path / "input.fa"
    write_fasta(fa, [(HDR_YAR, SEQ_A), (HDR_CAN, SEQ_B)])
    log.debug("input: 2 records with distinct sequences — no dups expected")
    run_dedup(tmp_path, [fa], extra_args=["--dedup_sequences"])
    exists = (tmp_path / "dedup_report.tsv").exists()
    log.debug("dedup_report.tsv present: %s (expected: False)", exists)
    assert not exists
    log.info("PASS: dedup_report.tsv not written when no seq dups")


# ---------------------------------------------------------------------------
# Multi-file merge preserves order and deduplicates across files
# ---------------------------------------------------------------------------

def test_cross_file_header_dedup(tmp_path):
    """
    A header appearing in file 1 and again in file 2 → only file 1 occurrence kept.
    Mirrors the pipeline merging relatives.fa + singletons.fa per query.
    """
    fa1 = tmp_path / "f1.fa"
    fa2 = tmp_path / "f2.fa"
    write_fasta(fa1, [(HDR_YAR, SEQ_A), (HDR_CAN, SEQ_B)])
    write_fasta(fa2, [(HDR_YAR, SEQ_C), (HDR_SCH, SEQ_C)])  # HDR_YAR is cross-file dup
    log.debug("f1: [YAR(SEQ_A), CAN(SEQ_B)]  f2: [YAR(SEQ_C), SCH(SEQ_C)] — YAR is cross-file dup")

    out, stdout = run_dedup(tmp_path, [fa1, fa2])
    headers = read_fasta_headers(out)
    log.debug("output headers (%d): %s | stdout: %s", len(headers), headers, stdout.strip())
    assert headers == [HDR_YAR, HDR_CAN, HDR_SCH]
    seqs = read_fasta_sequences(out)
    assert seqs[HDR_YAR] == SEQ_A   # file 1 version kept
    assert "header_dups_removed: 1" in stdout
    log.info("PASS: cross-file header dup removed — f1 occurrence of YAR kept")


def test_cross_file_seq_dedup(tmp_path):
    """
    Sequence appearing in file 1 under one header and file 2 under another →
    file 2 occurrence dropped when --dedup_sequences.
    """
    fa1 = tmp_path / "relatives.fa"
    fa2 = tmp_path / "extra.fa"
    write_fasta(fa1, [(HDR_YAR, SEQ_A), (HDR_CAN, SEQ_B)])
    write_fasta(fa2, [(HDR_ASP, SEQ_A), (HDR_NEU, SEQ_C)])  # HDR_ASP has same seq as HDR_YAR
    log.debug("relatives: [YAR(SEQ_A), CAN(SEQ_B)]  extra: [ASP(SEQ_A), NEU(SEQ_C)]")
    log.debug("ASP shares SEQ_A with YAR — should be dropped with --dedup_sequences")

    out, stdout = run_dedup(tmp_path, [fa1, fa2], extra_args=["--dedup_sequences"])
    headers = read_fasta_headers(out)
    log.debug("output headers (%d): %s | stdout: %s", len(headers), headers, stdout.strip())
    assert HDR_YAR in headers
    assert HDR_ASP not in headers   # dropped — same seq as HDR_YAR from file 1
    assert HDR_NEU in headers
    assert "seq_dups_removed: 1" in stdout
    log.info("PASS: cross-file seq dup removed — f1 YAR kept, f2 ASP dropped")


# ---------------------------------------------------------------------------
# Dedup with real pipeline sequences and deliberate duplicates
# ---------------------------------------------------------------------------

def test_real_sequences_duplicate_header_removed(tmp_path):
    """
    Duplicate header with real pipeline sequences: second occurrence dropped,
    first sequence preserved.
    yar02 appears twice with different sequences — only the first is kept.
    """
    fa = tmp_path / "input.fa"
    write_fasta(fa, [
        (HDR_YAR, REAL_SEQ_YAR02),
        (HDR_CAN, REAL_SEQ_CAN02),
        (HDR_YAR, REAL_SEQ_SCH02),   # duplicate header — different sequence
        (HDR_ASP, REAL_SEQ_ASP02),
    ])
    log.debug("input: 4 records, HDR_YAR appears twice with different real sequences")
    out, stdout = run_dedup(tmp_path, [fa])
    headers = read_fasta_headers(out)
    seqs = read_fasta_sequences(out)
    log.debug("output headers (%d): %s | stdout: %s", len(headers), headers, stdout.strip())
    assert len(headers) == 3
    assert headers.count(HDR_YAR) == 1
    assert seqs[HDR_YAR] == REAL_SEQ_YAR02   # first occurrence kept
    assert "header_dups_removed: 1" in stdout
    assert "seq_dups_removed: 0" in stdout
    log.info("PASS: real seq header dup removed — first occurrence of YAR kept")


def test_real_sequences_duplicate_sequence_removed(tmp_path):
    """
    Two proteins with different headers but identical sequences:
    second dropped only when --dedup_sequences is passed.
    Mirrors a real scenario where two species have the exact same protein.
    """
    fa = tmp_path / "input.fa"
    write_fasta(fa, [
        (HDR_YAR, REAL_SEQ_YAR02),
        (HDR_CAN, REAL_SEQ_CAN02),
        (HDR_SCH, REAL_SEQ_YAR02),   # different header, same seq as HDR_YAR
        (HDR_ASP, REAL_SEQ_ASP02),
        (HDR_NEU, REAL_SEQ_NEU02),
    ])
    log.debug("input: 5 records; SCH shares REAL_SEQ_YAR02 with YAR")

    # Without --dedup_sequences: all 5 retained
    out_no_dedup, _ = run_dedup(tmp_path, [fa])
    n_no_dedup = len(read_fasta_headers(out_no_dedup))
    log.debug("without --dedup_sequences: %d sequences retained (expected 5)", n_no_dedup)
    assert n_no_dedup == 5

    # With --dedup_sequences: HDR_SCH dropped (same seq as HDR_YAR)
    out_dedup, stdout = run_dedup(tmp_path, [fa], extra_args=["--dedup_sequences"])
    headers = read_fasta_headers(out_dedup)
    log.debug("with --dedup_sequences: %d sequences retained | stdout: %s", len(headers), stdout.strip())
    assert len(headers) == 4
    assert HDR_YAR in headers
    assert HDR_SCH not in headers
    assert "seq_dups_removed: 1" in stdout
    log.info("PASS: real seq dup removed with --dedup_sequences; YAR kept, SCH dropped")


def test_real_sequences_dedup_report_content(tmp_path):
    """
    dedup_report.tsv records the correct dropped/kept headers and sequence length
    when real pipeline sequences are used.
    """
    import pandas as pd
    fa = tmp_path / "input.fa"
    write_fasta(fa, [
        (HDR_YAR, REAL_SEQ_YAR02),
        (HDR_CAN, REAL_SEQ_YAR02),   # exact duplicate of REAL_SEQ_YAR02
        (HDR_SCH, REAL_SEQ_SCH02),
    ])
    log.debug("input: CAN shares REAL_SEQ_YAR02 with YAR (len=%d)", len(REAL_SEQ_YAR02))
    run_dedup(tmp_path, [fa], extra_args=["--dedup_sequences"])
    report = tmp_path / "dedup_report.tsv"
    assert report.exists()
    df = pd.read_csv(report, sep="\t")
    log.debug("report: %d rows | dropped=%s kept=%s len=%s",
              len(df), df.iloc[0]["dropped_header"], df.iloc[0]["kept_header"],
              df.iloc[0]["sequence_length"])
    assert len(df) == 1
    assert df.iloc[0]["dropped_header"] == HDR_CAN
    assert df.iloc[0]["kept_header"]    == HDR_YAR
    assert df.iloc[0]["sequence_length"] == len(REAL_SEQ_YAR02)
    log.info("PASS: report correctly records dropped=CAN kept=YAR len=%d", len(REAL_SEQ_YAR02))


def test_real_sequences_mixed_header_and_seq_dups(tmp_path):
    """
    Combined scenario: one header duplicate + one sequence duplicate.
    Both are removed independently.
    """
    fa = tmp_path / "input.fa"
    write_fasta(fa, [
        (HDR_YAR, REAL_SEQ_YAR02),
        (HDR_CAN, REAL_SEQ_CAN02),
        (HDR_YAR, REAL_SEQ_SCH02),   # header dup
        (HDR_NEU, REAL_SEQ_CAN02),   # seq dup (same as HDR_CAN)
        (HDR_ASP, REAL_SEQ_ASP02),
    ])
    log.debug("input: YAR appears twice (header dup); NEU shares seq with CAN (seq dup)")
    out, stdout = run_dedup(tmp_path, [fa], extra_args=["--dedup_sequences"])
    headers = read_fasta_headers(out)
    log.debug("output headers (%d): %s | stdout: %s", len(headers), headers, stdout.strip())
    assert len(headers) == 3   # YAR, CAN, ASP retained; SCH (header dup) and NEU (seq dup) dropped
    assert HDR_YAR in headers
    assert HDR_CAN in headers
    assert HDR_ASP in headers
    assert HDR_SCH not in headers   # header dup removed first
    assert HDR_NEU not in headers   # seq dup removed
    assert "header_dups_removed: 1" in stdout
    assert "seq_dups_removed: 1" in stdout
    log.info("PASS: mixed dups — header_dups=1 seq_dups=1; retained YAR+CAN+ASP")


# ---------------------------------------------------------------------------
# Dedup on real pipeline output (requires nextflow_test_dedup/ to exist)
#
# Run the pipeline first:
#   nextflow run workflows/echo.nf -params-file test/params.test.dedup.yaml
# Then run this file:
#   pytest test/test_dedup.py -v
# ---------------------------------------------------------------------------

import os as _os
import logging as _logging

_log       = _logging.getLogger("echo.test.dedup.pipeline")
_NF_DIR    = _os.path.abspath(_os.path.join(_os.path.dirname(__file__), ".."))
_DEDUP_DIR = _os.path.join(_NF_DIR, "nextflow_test_dedup")

_QUERY_SPECIES = ["saccharomyces_cerevisiae", "candida_albicans"]

# Expected sequence counts after the pipeline dedup run.
# Input data has 2 deliberate sequence duplicates (neu03≡sch03, yar02≡asp02) —
# both pairs appear in every query's relatives, so the pipeline removes 2 per query.
# Base count is 30 (29 ranked_cluster + 1 singleton); 30 - 2 = 28 after dedup.
_EXPECTED_SEQS = 28


@pytest.fixture(scope="module")
def dedup_outdir():
    if not _os.path.isdir(_DEDUP_DIR):
        pytest.skip(
            "nextflow_test_dedup/ not found — run the pipeline first:\n"
            "  nextflow run workflows/echo.nf -params-file test/params.test.dedup.yaml"
        )
    _log.info("dedup outdir: %s", _DEDUP_DIR)
    return _DEDUP_DIR


@pytest.mark.parametrize("query_name", _QUERY_SPECIES)
def test_pipeline_dedup_fasta_exists(dedup_outdir, query_name):
    """Pipeline output FASTA exists for each query."""
    path = _os.path.join(dedup_outdir, f"{query_name}_all_relatives.fa")
    _log.debug("checking FASTA exists: %s", path)
    assert _os.path.isfile(path), f"Missing: {path}"
    _log.info("PASS: %s_all_relatives.fa found", query_name)


@pytest.mark.parametrize("query_name", _QUERY_SPECIES)
def test_pipeline_dedup_sequence_count(dedup_outdir, query_name):
    """
    Pipeline removed the 2 deliberate sequence duplicates, leaving _EXPECTED_SEQS.
    Re-running dedup on the output finds zero further dups.
    """
    in_fa = _os.path.join(dedup_outdir, f"{query_name}_all_relatives.fa")
    headers = read_fasta_headers(in_fa)
    _log.debug("%s: found %d sequences in pipeline output (expected %d)",
               query_name, len(headers), _EXPECTED_SEQS)
    assert len(headers) == _EXPECTED_SEQS, \
        f"{query_name}: expected {_EXPECTED_SEQS} sequences, got {len(headers)}"
    _log.info("PASS: %s has %d sequences after pipeline dedup", query_name, len(headers))


@pytest.mark.parametrize("query_name", _QUERY_SPECIES)
def test_pipeline_dedup_output_is_clean(dedup_outdir, query_name, tmp_path):
    """Re-running dedup on the pipeline output finds no further duplicates."""
    in_fa = _os.path.join(dedup_outdir, f"{query_name}_all_relatives.fa")
    out_fa = tmp_path / f"{query_name}_rededup.fa"
    _log.debug("re-running dedup on %s -> %s", in_fa, out_fa)
    result = subprocess.run(
        [sys.executable, SCRIPT, "--inputs", in_fa, "--output", str(out_fa), "--dedup_sequences"],
        capture_output=True, text=True, cwd=str(tmp_path), check=False,
    )
    _log.debug("re-dedup stdout: %s", result.stdout.strip())
    assert result.returncode == 0, result.stderr
    assert "seq_dups_removed: 0" in result.stdout
    assert "header_dups_removed: 0" in result.stdout
    assert not (tmp_path / "dedup_report.tsv").exists()
    _log.info("PASS: %s pipeline output is clean — no residual dups", query_name)


def test_pipeline_dedup_report_exists(dedup_outdir):
    """dedup_report.tsv is published to the output dir when seq dups exist."""
    report = _os.path.join(dedup_outdir, "dedup_report.tsv")
    _log.debug("checking dedup_report.tsv exists: %s", report)
    assert _os.path.isfile(report), "dedup_report.tsv missing — pipeline did not detect seq dups"
    _log.info("PASS: dedup_report.tsv found in pipeline output dir")


def test_pipeline_dedup_report_content(dedup_outdir):
    """
    dedup_report.tsv has exactly 2 rows — one per deliberate duplicate pair.
    Pair 1: neu03 / sch03  (sequence_length = 48)
    Pair 2: yar02 / asp02  (sequence_length = 60)
    """
    import pandas as pd
    report = _os.path.join(dedup_outdir, "dedup_report.tsv")
    df = pd.read_csv(report, sep="\t")
    _log.debug("report: %d rows, columns: %s", len(df), list(df.columns))
    _log.debug("report contents:\n%s", df.to_string(index=False))
    assert list(df.columns) == ["dropped_header", "kept_header", "sequence_length"]
    assert len(df) == 2, f"Expected 2 rows in dedup_report.tsv, got {len(df)}"

    lengths = sorted(df["sequence_length"].tolist())
    _log.debug("dup sequence lengths (sorted): %s", lengths)
    assert lengths == [48, 60], f"Unexpected sequence lengths: {lengths}"

    # Pair 1: neu03 / sch03 (len 48)
    row48 = df[df["sequence_length"] == 48].iloc[0]
    pair48 = {row48["dropped_header"].split("|")[0], row48["kept_header"].split("|")[0]}
    _log.debug("len-48 pair: %s  (dropped=%s kept=%s)",
               pair48, row48["dropped_header"], row48["kept_header"])
    assert pair48 == {"neu03", "sch03"}, f"Unexpected pair for len-48 dup: {pair48}"

    # Pair 2: yar02 / asp02 (len 60)
    row60 = df[df["sequence_length"] == 60].iloc[0]
    pair60 = {row60["dropped_header"].split("|")[0], row60["kept_header"].split("|")[0]}
    _log.debug("len-60 pair: %s  (dropped=%s kept=%s)",
               pair60, row60["dropped_header"], row60["kept_header"])
    assert pair60 == {"yar02", "asp02"}, f"Unexpected pair for len-60 dup: {pair60}"
    _log.info("PASS: dedup_report.tsv has correct 2 pairs — (neu03,sch03) and (yar02,asp02)")
