# Pomegranate Description Plot Builder

This repository contains an R script for creating pomegranate annotation description files used in description plots.

The script builds a merged annotation table from UniProt data, RefSeq identifiers, and Gene Ontology (GO) information. The resulting file is intended for downstream visualization and functional summary workflows related to pomegranate datasets.

## Purpose

The goal of this repository is to prepare a clean description file for pomegranate genes or proteins by combining several annotation sources into one output table.

## Features

- Extracts protein descriptions from UniProt data
- Merges UniProt annotations with RefSeq IDs
- Adds Gene Ontology information
- Produces a description file suitable for pomegranate description plots
- Helps standardize annotation input for downstream analysis

## Input

The script is designed to work with annotation files containing:

- UniProt protein annotation data
- RefSeq ID mappings
- GO annotation data

## Output

The output is a merged description file that may include:

- UniProt accession IDs
- Protein names or descriptions
- RefSeq identifiers
- Gene names, when available
- GO terms or GO annotations

## Usage

Run the R script:

```r
source("build_uniprot_description_file.r")
