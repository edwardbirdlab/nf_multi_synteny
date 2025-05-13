#!/bin/bash

# Usage: ./replace_ids.sh SequenceIDs.txt input.blast output.blast

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 SequenceIDs.txt input.blast output.blast"
    exit 1
fi

SEQIDS="$1"
BLAST="$2"
OUT="$3"

awk '
    # Read mapping file
    FNR==NR {
        split($0, a, ":")
        id = a[1]
        split(a[2], b, " ")
        transcript = b[1]
        map[id] = transcript
        next
    }

    # Replace columns 1 and 2
    {
        if ($1 in map) $1 = map[$1]
        if ($2 in map) $2 = map[$2]
        print
    }
' "$SEQIDS" "$BLAST" > "$OUT"

echo "✔️ Replaced IDs written to $OUT"
