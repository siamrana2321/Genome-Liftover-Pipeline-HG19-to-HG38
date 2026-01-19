#!/usr/bin/env python3
"""
validate.py

Validation for liftover_project outputs:
- checks required columns
- checks ncbi_build == GRCh38
- checks REF allele matches FASTA at Chrom:Start
- counts SNV vs INDEL mismatches
- writes a summary JSON and optional CSV for details
"""

import csv
import json
from pathlib import Path
from collections import Counter
from pyfaidx import Fasta

# CONFIG: adjust paths / thresholds here or read from config.yaml if you prefer
PROJECT_ROOT = Path(__file__).resolve().parent
RESOURCES = PROJECT_ROOT / "resources"
DATA_OUT = PROJECT_ROOT / "data" / "output"
UNMAP_DIR = PROJECT_ROOT / "data" / "unmap"
SUMMARY_DIR = PROJECT_ROOT / "logs"
FASTA = RESOURCES / "GRCh38.primary_assembly.genome.fa"

# validation thresholds
MAX_REF_MISMATCH_RATE = 0.01  # 1% mismatches allowed before failing
FAIL_ON_HIGH_MISMATCH = False  # set True to exit non-zero if above threshold

# MAF/CBioPortal columns we expect (common names)
REQ_COLS = [
    "Chromosome",
    "Start_Position",
    "End_Position",
    "Reference_Allele",
    "Tumor_Seq_Allele2",
    "NCBI_Build"
]


def read_header_and_rows(file_path):
    """Return (header_list, rows_list). Skips initial comment lines starting with #."""
    with open(file_path, "r", newline="") as fh:
        all_lines = fh.readlines()

    # find first non-comment line (the header)
    start_idx = 0
    for i, line in enumerate(all_lines):
        if not line.startswith("#"):
            start_idx = i
            break

    # create a csv.DictReader over the header+rows portion
    reader = csv.DictReader(all_lines[start_idx:], delimiter="\t")
    header = reader.fieldnames
    rows = list(reader)  # materialize to avoid closed-file problems
    return header, rows


def check_columns(header):
    missing = [c for c in REQ_COLS if c not in header]
    return missing


def validate_file(file_path, fasta_obj):
    header, rows = read_header_and_rows(file_path)
    missing = check_columns(header)
    if missing:
        return {
            "file": str(file_path),
            "error": "missing_columns",
            "missing_columns": missing
        }

    counts = Counter()
    bad_records = []
    total = 0
    for row in rows:
        total += 1
        try:
            chrom = row["Chromosome"].strip()
            start = int(row["Start_Position"])
            ref = row["Reference_Allele"].strip()
            alt = row["Tumor_Seq_Allele2"].strip()
            build = row.get("ncbi_build", "").strip()
        except Exception as e:
            counts["parse_error"] += 1
            bad_records.append({"reason": "parse_error", "error": str(e), "row": row})
            continue

        counts["total"] += 1
        if build and "GRCh38" not in build:
            counts["wrong_build"] += 1

        # normalize chromosome name for FASTA lookup if necessary
        seqname = chrom
        if seqname not in fasta_obj:
            # try with chr prefix, or without
            if chrom.startswith("chr"):
                seqname = chrom[3:]
            else:
                seqname = "chr" + chrom
        if seqname not in fasta_obj:
            counts["no_chrom_in_fasta"] += 1
            bad_records.append({
                "reason": "chrom_missing_in_fasta",
                "chrom": chrom,
                "start": start,
                "ref": ref,
                "alt": alt
            })
            continue

        # fetch reference base(s) from fasta
        try:
            # MAF positions are 1-based inclusive; pyfaidx uses 1-based inclusive via slice indices
            fetched = fasta_obj[seqname][start - 1 : start - 1 + len(ref)].seq
        except Exception as e:
            counts["fasta_fetch_error"] += 1
            bad_records.append({
                "reason": "fasta_fetch_error",
                "chrom": seqname,
                "start": start,
                "error": str(e)
            })
            continue

        if fetched.upper() != ref.upper():
            counts["ref_mismatch"] += 1
            # differentiate SNV vs indel
            if len(ref) == 1 and len(alt) == 1:
                counts["snv_ref_mismatch"] += 1
            else:
                counts["indel_ref_mismatch"] += 1

            bad_records.append({
                "reason": "ref_mismatch",
                "chrom": seqname,
                "start": start,
                "ref_expected": fetched,
                "ref_observed": ref,
                "alt": alt
            })
        else:
            counts["ref_ok"] += 1

    mismatch_rate = counts.get("ref_mismatch", 0) / max(counts.get("total", 0), 1)

    result = {
        "file": str(file_path),
        "total_variants": counts.get("total", 0),
        "ref_ok": counts.get("ref_ok", 0),
        "ref_mismatch": counts.get("ref_mismatch", 0),
        "snv_ref_mismatch": counts.get("snv_ref_mismatch", 0),
        "indel_ref_mismatch": counts.get("indel_ref_mismatch", 0),
        "wrong_build": counts.get("wrong_build", 0),
        "no_chrom_in_fasta": counts.get("no_chrom_in_fasta", 0),
        "mismatch_rate": mismatch_rate,
        "bad_records_sample": bad_records[:50]
    }
    return result


def main():
    fasta = Fasta(str(FASTA), sequence_always_upper=True)
    SUMMARY_DIR.mkdir(parents=True, exist_ok=True)
    reports = []
    for f in sorted(DATA_OUT.glob("*.txt")):
        print("Validating:", f.name)
        report = validate_file(f, fasta)
        reports.append(report)
        out_json = SUMMARY_DIR / f"{f.stem}.validation.json"
        out_json.write_text(json.dumps(report, indent=2))

    agg = {
        "files_validated": len(reports),
        "reports": reports
    }
    (SUMMARY_DIR / "validation_summary.json").write_text(json.dumps(agg, indent=2))
    print("Validation complete. Summary written to", SUMMARY_DIR / "validation_summary.json")

    overall_mismatch = max((r.get("mismatch_rate", 0) for r in reports), default=0)
    if FAIL_ON_HIGH_MISMATCH and overall_mismatch > MAX_REF_MISMATCH_RATE:
        print("ERROR: mismatch rate above threshold:", overall_mismatch)
        raise SystemExit(2)


if __name__ == "__main__":
    main()
