import subprocess
import shutil
import yaml
from pathlib import Path
import pandas as pd
from io import StringIO

# ================= PATH SETUP =================
BASE_DIR = Path(__file__).resolve().parent
RESOURCES = BASE_DIR / "resources"
DATA = BASE_DIR / "data"
INPUT = DATA / "input"
OUTPUT = DATA / "output"
UNMAP = DATA / "unmap"
LOGS = BASE_DIR / "logs"
CONFIG = BASE_DIR / "config.yaml"
# ==============================================

with open(CONFIG) as f:
    cfg = yaml.safe_load(f)

CHAIN = RESOURCES / "hg19ToHg38.over.chain.gz"
FASTA = RESOURCES / "GRCh38.primary_assembly.genome.fa"
BUILD = cfg["reference"]["build"]
CHROM_STYLE = cfg["chromosomes"]["style"]
WIN_OUTPUT = Path(cfg["windows"]["output_dir"])

# ================= FINAL OUTPUT SCHEMA =================
FINAL_COLUMNS = [
    "Hugo_Symbol",
    "Entrez_Gene_Id",
    "NCBI_Build",
    "Chromosome",
    "Start_Position",
    "End_Position",
    "Consequence",
    "Variant_Classification",
    "Variant_Type",
    "Reference_Allele",
    "Tumor_Seq_Allele1",
    "Tumor_Seq_Allele2",
    "Tumor_Sample_Barcode",
    "Transcript_ID",
    "RefSeq",
    "Gene",
    "Annotation_Status",
    "Filter",
    "Tissue",
    "Cancer_Type",
    "PMID",
    "Study",
    "Seq_Tech"
]

# ======================================================
# SAFE Variant_Type recomputation + column enforcement
# ======================================================
def process_output_maf(maf_path: Path):
    # Read raw file
    with open(maf_path, "r") as f:
        lines = f.readlines()

    # Preserve CrossMap comments
    comment_lines = [l for l in lines if l.startswith("#")]
    table_lines = [l for l in lines if not l.startswith("#")]

    if not table_lines:
        raise RuntimeError("No data found in MAF file")

    df = pd.read_csv(
        StringIO("".join(table_lines)),
        sep="\t",
        dtype=str,
        keep_default_na=False,
        na_filter=False
    )

    # -------- Recompute Variant_Type --------
    def classify(row):
        ref = row.get("Reference_Allele", "")
        alt = row.get("Tumor_Seq_Allele2", "")
        if ref == "-" and alt == "-":
            return "-"
        if alt == "-" and ref != "-":
            return "DEL"
        if ref == "-" and alt != "-":
            return "INS"
        if len(ref) > len(alt):
            return "DEL"
        if len(ref) < len(alt):
            return "INS"
        if len(ref) == len(alt) == 1:
            return "SNP"
        if len(ref) == len(alt) > 1:
            return "ONP"
        return "-"

    df["Variant_Type"] = df.apply(classify, axis=1)

    # -------- Enforce final schema --------
    for col in FINAL_COLUMNS:
        if col not in df.columns:
            df[col] = "-"

    df = df[FINAL_COLUMNS]
    df.replace("", "-", inplace=True)

    # -------- Write back TXT safely --------
    with open(maf_path, "w") as f:
        for l in comment_lines:
            f.write(l)
        df.to_csv(f, sep="\t", index=False)

    return df


# ======================================================
# CrossMap runner
# ======================================================
def run_crossmap(input_file: Path):
    OUTPUT.mkdir(exist_ok=True)
    UNMAP.mkdir(exist_ok=True)
    LOGS.mkdir(exist_ok=True)
    WIN_OUTPUT.mkdir(parents=True, exist_ok=True)

    output_file = OUTPUT / f"{input_file.stem}.{BUILD}.txt"
    log_file = LOGS / f"{input_file.stem}.log"

    cmd = [
        "CrossMap", "maf",
        "--chromid", CHROM_STYLE,
        str(CHAIN),
        str(input_file),
        str(FASTA),
        BUILD,
        str(output_file)
    ]

    print(f"\n▶ Liftover started: {input_file.name}")

    with open(log_file, "w") as log:
        subprocess.run(cmd, stdout=log, stderr=log, check=True)

    # Post-process (Variant_Type + column filtering)
    df = process_output_maf(output_file)

    # Handle unmapped variants
    unmap_file = Path(str(output_file) + ".unmap")
    if unmap_file.exists():
        shutil.move(unmap_file, UNMAP / unmap_file.name)

    # Copy TXT to Windows
    win_txt = WIN_OUTPUT / output_file.name
    shutil.copy(output_file, win_txt)

    # Write CSV to Windows
    win_csv = WIN_OUTPUT / f"{output_file.stem}.csv"
    df.to_csv(win_csv, index=False)

    print(f"✔ Completed: {output_file.name}")
    print(f"✔ Windows TXT: {win_txt.name}")
    print(f"✔ Windows CSV: {win_csv.name}")


# ======================================================
# Main
# ======================================================
def main():
    INPUT.mkdir(parents=True, exist_ok=True)
    files = list(INPUT.glob("*.txt"))

    if not files:
        print("❌ No input files found in data/input/")
        return

    for f in files:
        run_crossmap(f)


if __name__ == "__main__":
    main()
