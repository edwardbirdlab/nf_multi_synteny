#!/bin/bash

# Usage: ./rename_fasta_headers.sh mapping.tsv input.fasta output.fasta

MAPFILE=$1
INPUT_FASTA=$2
OUTPUT_FASTA=$3

if [[ -z "$MAPFILE" || -z "$INPUT_FASTA" || -z "$OUTPUT_FASTA" ]]; then
    echo "Usage: $0 mapping.tsv input.fasta output.fasta"
    echo "mapping.tsv should be a 2-column tab-separated file: new_name<TAB>old_name"
    exit 1
fi

awk 'NR==FNR {
    map[$2]=$1;
    next
}
/^>/ {
    header=$0;
    sub(/^>/,"",header);
    if (header in map)
        print ">"map[header];
    else
        print ">"header;
    next
}
{ print }
' "$MAPFILE" "$INPUT_FASTA" > "$OUTPUT_FASTA"

echo "Renamed FASTA saved to: $OUTPUT_FASTA"
