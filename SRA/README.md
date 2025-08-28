# SRA ONT Open Data Fetching Pipeline

A simple bioinformatics workflow designed for automated retrieval, filtration, and annotation of publicly available Oxford Nanopore Technologies (ONT) sequencing data from NCBI's Sequence Read Archive (SRA).

## ✨ Features

- **Data Fetching**: Automated retrieval of ONT sequencing data from SRA public repositories
- **Multi-level Filtration**: Customizable filtering pipeline to select datasets meeting specific research criteria
- **Comprehensive Annotation**: Enrichment with essential metadata including:
  - Raw sequencing data formats
  - BioSample attributes
  - Reference genome sequences for corresponding species
  - Taxonomic classification and phylogenetic context

- **This pipeline is mainly based on Entrez Direct (EDirect) by recruiting different parts of the information and filtering them.**

## 🚀 Quick Start

### Prerequisites
- Conda package manager (Miniconda or Anaconda)
- Internet connectivity for SRA database access

### Installation & Execution

1. **Create the computational environment**:
```bash
conda env create -f environment.yml
```
2. **Activate the environment**:
```bash
conda activate ont_datafetch
```
3. **Execute the pipeline**:
```bash
python3 run_datafetch.py -o OUTPUT_DIR [-f]
```
**NOTICE** that the progress could run hours to days

4. **Check the results:
```bash
sra_runinfo.txt                             # Raw runinfo get from esearch and efetch
sra_runinfo_filtered.tsv                    # Filtered runinfo
sra_runinfo_filtered_with_rawformat.tsv     # Add raw uploaded data format to the table
assembly_summary_genbank.txt                # GenBank assembly summary file
sra_runinfo_integrated.tsv                  # Integrating the runs by biosample and assign a refseq
```
**The results need more mannually checking and filtering.**
