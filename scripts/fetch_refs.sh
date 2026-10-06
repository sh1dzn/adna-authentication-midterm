#!/usr/bin/env bash
# Fetch the two reference sequences from NCBI. Accessions are the provenance.
#   NC_003143.1  Yersinia pestis CO92 chromosome (4,653,728 bp)
#   NC_012920.1  Homo sapiens mitochondrion, rCRS (16,569 bp)
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p data/ref
E="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&retmode=text&id"
[ -s data/ref/ypestis_co92.fa ] || curl -fsS --retry 4 "$E=NC_003143.1" -o data/ref/ypestis_co92.fa
[ -s data/ref/human_rCRS.fa ]   || curl -fsS --retry 4 "$E=NC_012920.1" -o data/ref/human_rCRS.fa
grep -c '>' data/ref/*.fa
