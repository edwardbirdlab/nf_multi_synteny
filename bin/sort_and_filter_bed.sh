#!/bin/bash

# Usage: ./sort_and_filter_bed.sh combined.bed output.bed

# Check for input and output arguments
if [[ $# -ne 2 ]]; then
    echo "Usage: $0 input.bed output.bed"
    exit 1
fi

INPUT_BED=$1
OUTPUT_BED=$2

# Check if the input file exists
if [[ ! -f "$INPUT_BED" ]]; then
    echo "Error: File $INPUT_BED not found!"
    exit 1
fi

# Sort by chromosome and start position, then reorder columns to 1,4,2,3
sort -k1,1 -k2,2n "$INPUT_BED" | awk 'BEGIN {OFS="\t"} $4 != "" {print $1, $4, $2, $3}' > "$OUTPUT_BED"

# Check if the output file was created
if [[ -f "$OUTPUT_BED" ]]; then
    echo "Sorting and filtering completed. Output saved to $OUTPUT_BED."
else
    echo "Error: Output file not created."
fi
