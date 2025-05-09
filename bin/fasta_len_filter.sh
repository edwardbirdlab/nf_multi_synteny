#!/bin/bash

# Usage: ./filter_fasta_by_length.sh input.fasta MIN_LENGTH output.fasta

INPUT=$1
MINLEN=$2
OUTPUT=$3

if [[ -z "$INPUT" || -z "$MINLEN" || -z "$OUTPUT" ]]; then
    echo "Usage: $0 input.fasta MIN_LENGTH output.fasta"
    exit 1
fi

awk -v minlen="$MINLEN" '
BEGIN { RS=">"; ORS=""; }
NR > 1 {
    header = $1;
    sub(/[^\n]*\n/, "", $0);  # remove header line
    gsub(/\n/, "", $0);       # join all sequence lines
    if (length($0) >= minlen) {
        print ">" header "\n";
        for (i = 1; i <= length($0); i += 60)
            print substr($0, i, 60) "\n";
    }
}
' "$INPUT" > "$OUTPUT"

echo "Filtered FASTA saved to: $OUTPUT"