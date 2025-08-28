# Find ONT Open Data

A collection of scripts for retrieving ONT whole-genome sequencing data from public databases for what we interested in.

## Overview

This repository contains practical scripts and workflows used to query, identify, and download Oxford Nanopore Technologies (ONT) whole-genome sequencing data from various public genomic databases. These scripts may help for efficiently access animal genomic data relevant to our interests.

## Purpose

The primary goal is to provide simple, reusable scripts that:
- Search public sequence databases (SRA, ENA, etc.) for ONT WGS data
- Filter results based on basic criteria relevant to our projects
- Download datasets(not sequences) and associated metadata
- Organize retrieved data in a consistent structure

## Repository Structure
```plaintext
./
├── genomeark/     # Scripts for querying and downloading from GenomeArk
├── GSA/           # Scripts for working with Genome Sequence Archive (GSA)
├── org.one/       # Scripts for org.one project
├── QC/            # Quality control using LongBow
├── SRA/           # Scripts for NCBI Sequence Read Archive access
└── README.md      # This documentation file
```

