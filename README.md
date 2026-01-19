# Genome Liftover Pipeline GRCh37 → GRCh38 (cBioPortal / MAF)

## 📌 Overview

This project implements a **fully automated, reproducible pipeline** to lift genomic mutation data from **GRCh37 (hg19)** to **GRCh38 (hg38)** using **CrossMap** within **WSL (Windows Subsystem for Linux)**.

The pipeline is designed for **cBioPortal / MAF-like mutation files** and ensures:

- Correct coordinate liftover (GRCh37 → GRCh38)
- Biological correctness through allele-aware processing
- Recalculation of allele-derived fields (e.g., `Variant_Type`)
- Standardization to a fixed, predefined schema
- Dual-format output for both computational and reporting use (`.txt` and `.csv`)
- Safe handling of missing annotations

This repository is suitable for **research, clinical bioinformatics, and production-grade data processing**.

---

## 🧬 Key Features

- ✅ Uses **CrossMap (maf mode)** for accurate liftover
- ✅ Validates and recomputes `Variant_Type` after liftover
- ✅ Preserves only required columns in final output
- ✅ Automatically adds missing columns with placeholder values (`-`)
- ✅ Outputs:
  - Tab-delimited `.txt` (for bioinformatics tools)
  - CSV `.csv` (for Excel, Power BI, reporting)
- ✅ Compatible with **cBioPortal** ingestion
- ✅ Designed to run on **Windows via WSL**

---

## 📂 Project Structure

```text
liftover_project/
├── liftover.py              # Main pipeline script
├── config.yaml              # Configuration file
├── resources/               # Reference files
│   ├── hg19ToHg38.over.chain.gz
│   └── GRCh38.primary_assembly.genome.fa
├── data/
│   ├── input/               # Input GRCh37 MAF-like .txt files
│   ├── output/              # Internal GRCh38 outputs
│   └── unmap/               # Unmapped variants
├── logs/                    # CrossMap logs
└── README.md

## 🧾 Final Output Schema
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

Note :
If a column exists → its data is preserved
If a column is missing or empty → it is created and filled with "-"

## 🛠️ Requirements
Operating System

Windows 10 / 11 with WSL2
Software (inside WSL)
Ubuntu (recommended)
Python ≥ 3.10
CrossMap
pandas
PyYAML
