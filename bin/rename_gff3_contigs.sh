#!/bin/bash

# Usage: ./rename_gff3_contigs.sh mapping.tsv input.gff3 output.gff3
# mapping.tsv should be a 2-column tab-separated file: new_name<TAB>old_full_header

MAPFILE=$1
INPUT_GFF=$2
OUTPUT_GFF=$3

if [[ -z "$MAPFILE" || -z "$INPUT_GFF" || -z "$OUTPUT_GFF" ]]; then
    echo "Usage: $0 mapping.tsv input.gff3 output.gff3"
    echo "mapping.tsv should be a 2-column tab-separated file: new_name<TAB>old_full_header"
    exit 1
fi

awk 'BEGIN { FS=OFS="\t" }
NR==FNR {
    split($2, a, " ");
    map[a[1]] = $1;
    next
}
/^#/ {
    print;
    next
}
($1 in map) {
    $1 = map[$1];
    print;
}
' "$MAPFILE" "$INPUT_GFF" > "$OUTPUT_GFF"

echo "Renamed + filtered GFF3 saved to: $OUTPUT_GFF"
