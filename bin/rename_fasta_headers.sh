#!/bin/bash

# Usage: ./rename_fasta_headers.sh mapping.tsv input.fasta output.fasta
# mapping.tsv should be a 2-column tab-separated file: new_name<TAB>old_short_name

MAPFILE=$1
INPUT_FASTA=$2
OUTPUT_FASTA=$3

if [[ -z "$MAPFILE" || -z "$INPUT_FASTA" || -z "$OUTPUT_FASTA" ]]; then
    echo "Usage: $0 mapping.tsv input.fasta output.fasta"
    echo "mapping.tsv should be a 2-column tab-separated file: new_name<TAB>old_short_name"
    exit 1
fi

awk 'NR==FNR {
    map[$2]=$1;
    next
}
/^>/ {
    header=$0;
    sub(/^>/,"",header);
    split(header, a, " ");
    key = a[1];
    if (key in map)
        print ">" map[key];
    else
        print ">" header;
    next
}
{ print }
' "$MAPFILE" "$INPUT_FASTA" > "$OUTPUT_FASTA"

echo "Renamed FASTA saved to: $OUTPUT_FASTA"
