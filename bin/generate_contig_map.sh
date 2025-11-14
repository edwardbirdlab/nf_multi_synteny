#!/bin/bash

# Usage: ./generate_contig_map.sh input.fasta output.tsv PREFIX

INPUT_FASTA=$1
OUTPUT_TSV=$2
PREFIX=$3

if [[ -z "$INPUT_FASTA" || -z "$OUTPUT_TSV" || -z "$PREFIX" ]]; then
    echo "Usage: $0 input.fasta output.tsv PREFIX"
    echo "Example: $0 genome.fasta genome_map.tsv ctg_"
    exit 1
fi

grep '^>' "$INPUT_FASTA" | cut -c2- | while read -r header; do
    full="$header"
    short=$(echo "$header" | awk '{print $1}')
    echo -e "${PREFIX}${short}\t${short}"
done > "$OUTPUT_TSV"

echo "Contig map saved to: $OUTPUT_TSV"
