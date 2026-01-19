# Genome Liftover Pipeline GRCh37 → GRCh38 (cBioPortal / MAF)

## Overview

This project implements a **fully automated, reproducible pipeline** to convert cancer mutation data from **GRCh37 (hg19)** to **GRCh38 (hg38)** using **CrossMap**.  
It is designed for **cBioPortal / MAF-like mutation files** and includes:

- Coordinate liftover (GRCh37 → GRCh38)
- Allele-aware validation using the GRCh38 reference genome
- Recalculation of allele-derived fields (e.g., `Variant_Type`)
- Schema standardization to a fixed set of columns
- Automatic handling of missing columns and values
- Dual-format output for Windows users: **TXT (tab-delimited)** and **CSV**
- Built-in validation to ensure biological correctness
- Clear logging and unmapped-variant tracking

The pipeline runs on **Windows using WSL (Ubuntu)** and is suitable for research, production, and audit-reviewed environments.

---

## Key Features

- Uses UCSC hg19→hg38 chain files via CrossMap
- Preserves biological correctness (allele-aware liftover)
- Recomputes `Variant_Type` according to MAF specification
- Enforces a fixed, documented output schema
- Automatically adds missing columns with placeholder values (`-`)
- Generates both `.txt` and `.csv` outputs for downstream use
- Includes post-liftover validation checks
- Compatible with cBioPortal and downstream annotation tools

---

## Final Output Schema

Both TXT and CSV outputs contain **only** the following columns (in this order):

Hugo_Symbol  
Entrez_Gene_Id  
NCBI_Build  
Chromosome  
Start_Position  
End_Position  
Consequence  
Variant_Classification  
Variant_Type  
Reference_Allele  
Tumor_Seq_Allele1  
Tumor_Seq_Allele2  
Tumor_Sample_Barcode  
Transcript_ID  
RefSeq  
Gene  
Annotation_Status  
Filter  
Tissue  
Cancer_Type  
PMID  
Study  
Seq_Tech  

Missing data or missing columns are filled with `-`.

---

## System Requirements

- Windows 10 / 11
- Windows Subsystem for Linux (WSL2 recommended)
- Ubuntu (via WSL)
- Python 3.10+
- CrossMap
- samtools
- pandas
- pyfaidx

---

## Installation (From Scratch)

### 1. Install WSL (PowerShell, Admin)
```powershell
wsl --install
```

---

### 2. Install Linux dependencies
```bash
sudo apt update && sudo apt upgrade -y
sudo apt install -y wget gzip build-essential python3-venv python3-pip samtools
```

---

### 3. Create Python environment
```bash
python3 -m venv crossmap_env
source crossmap_env/bin/activate
pip install --upgrade pip
pip install CrossMap pandas pyfaidx pyyaml
```

---

### 4. Create project structure
```bash
mkdir -p liftover_project/{resources,data/input,data/output,data/unmap,logs}
cd liftover_project
```

---

### 5. Download reference files

UCSC chain file (hg19 → hg38):
```bash
wget -P resources/ https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz
```

GRCh38 FASTA (GENCODE):
Download from https://www.gencodegenes.org/human/

Then:
```bash
samtools faidx resources/GRCh38.primary_assembly.genome.fa
```

---

### 6. Configure project (`config.yaml`)

```yaml
reference:
  build: GRCh38

chromosomes:
  style: s

windows:
  output_dir: /mnt/d/Filtered_CBiosportal_Output
```

---

## Running the Pipeline

### Copy input file
```bash
cp "/mnt/d/Path/To/Input/data_mutations.txt" liftover_project/data/input/
```

### Run liftover
```bash
cd liftover_project
source ~/crossmap_env/bin/activate
python liftover.py
python validate.py
```

---

## Validation and Quality Control

After liftover and post-processing, the project provides a dedicated validation step implemented in `validate.py` to verify the biological and structural correctness of GRCh38 outputs.

The validation script is designed to detect silent errors that may occur during coordinate liftover or allele normalization.

### Validation checks performed:

**1. Required column checks**

  The validator confirms the presence of essential MAF/cBioPortal fields, including:

  - Chromosome
  - Start_Position
  - End_Position
  - Reference_Allele
  - Tumor_Seq_Allele2
  - NCBI_Build

  Files missing required columns are flagged immediately.

**2. Reference build verification**

  Each record is checked to ensure that `NCBI_Build` corresponds to GRCh38.

**3. Reference allele verification against FASTA**

  For each variant:
  - The reference allele is fetched directly from the GRCh38 FASTA
  - The fetched sequence is compared against `Reference_Allele`
  - Mismatches are recorded and categorized

  This ensures that lifted coordinates remain biologically consistent with the target genome.

**4. SNV vs INDEL mismatch classification**

  Reference mismatches are further classified into:
  - SNV reference mismatches
  - INDEL reference mismatches

  This helps distinguish systematic coordinate errors from localized alignment issues.

**5. Chromosome availability checks**
  Variants mapping to chromosomes not present in the reference FASTA are detected and reported.

**6. Detailed reporting and audit trail**
  For each processed file, the validator generates:

  - A per-file JSON report in logs/
  - An aggregated validation_summary.json
  - A sample of problematic records for manual inspection

  Validation thresholds (e.g., maximum allowed mismatch rate) are configurable directly in validate.py.
---

## Output

### Inside WSL
- `data/output/` → GRCh38 tab-delimited files
- `data/unmap/` → unmapped variants
- `logs/` → execution logs

### Windows Output
```
*.GRCh38.txt
*.GRCh38.csv
```

---

## Notes on Variant_Type

`Variant_Type` is an **allele-derived field** and is not preserved blindly.  
After liftover, it is recomputed according to the MAF specification:

- `DEL` – deletion
- `INS` – insertion
- `SNP` – single-nucleotide polymorphism
- `ONP` – multi-nucleotide polymorphism

This ensures biological correctness and downstream compatibility.

---

## Best Practices

- Always inspect `data/unmap/` after a run
- Avoid opening raw MAF files in Excel for validation
- Use TSV-aware tools (`awk`, `pandas`, `cut`) for inspection
- Keep reference FASTA and chain versions documented

---

## License

This project is intended for academic and research use.
Please ensure compliance with data governance and privacy policies.

---

## Author

Mohammad Siam Ahmed Rana


