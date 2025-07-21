#!/usr/bin/env bash

count_files=("$1"/*/*gene_counts.txt)
read_counts="${count_files[0]}"
cut -f1 "$read_counts" > gene_read_counts.txt
for cf in "${count_files[@]}"; do
    paste gene_read_counts.txt <(cut -f2- "$cf") > merged.tmp && mv merged.tmp gene_read_counts.txt
done
header=$(printf "Gene_ID\t%s\n" "$(basename -a "${count_files[@]%%.*}")" | tr '\n' '\t' | sed 's/\t$/\n/')
{ echo "$header"; cat gene_read_counts.txt; } > merged.tmp && mv merged.tmp gene_read_counts.txt
