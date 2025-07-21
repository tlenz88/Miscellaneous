#!/bin/bash

# List of sample names
samples1=("HFF")
samples2=("ME49RFP_6hpi", "ME49RFP_24hpi")

# List of chromosome numbers
chromosomes=("01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12" "13" "14", "15", "16", "17", "18", "19", "20", "21", "22")
for i in "${samples1[@]}"; do
    for j in "${samples2[@]}"; do
        for k in "${chromosomes[@]}"; do
            python3 /mnt/f/scripts/hic_scripts/selfish.py \
                -m1 "${i}/${i}_100000_iced_cpm.matrix" \
                -m2 "${j}/${j}_100000_iced_cpm.matrix" \
                -ch "chr${k}" \
                -r 100000 \
                -bed1 "${i}/${i}_100000_abs.bed" \
                -bed2 "${j}/${j}_100000_abs.bed" \
                -o "${i}_vs_${j}_chr${k}.txt" \
                -t 0.05
        done
    done
done

echo "Merging differential interaction matrices."
prefixes=$(find $1 -type f -name '*.txt' -exec basename {} \; | awk -F_ '{print $1"_"$2}' | sort | uniq)
for prefix in ${prefixes[@]}; do
    output_file="${prefix}_diff.txt"
    echo "Saving ${prefix} matrices to ${output_file}"
    first_file=true
    for file in ${prefix}_*.txt; do
        if $first_file; then
            cat "$file" > "$output_file"
            first_file=false
        else
            tail -n +2 "$file" >> "$output_file"
        fi
    done
done

prefixes=$(find $1 -type f -name '*.txt' -exec basename {} \; | awk -F_ '{print $1"_"$2}' | sort | uniq)
for prefix in ${prefixes[@]}; do
    echo $prefix
    python3 /mnt/f/scripts/hic_scripts/plot_selfish.py $prefix /mnt/f/organism_genome/CparvumIowaII/CryptoDB-68_CparvumIowaII_Genome.chrom.sizes ${prefix}_diff.txt
done
: << EOF
#prefixes=$(find $1 -type f -name '*.txt' -exec basename {} \; | awk -F_ '{print $1"_"$2}' | sort | uniq)
for prefix in ${prefixes[@]}; do
    echo "Plotting gene differential interaction plot for ${prefix}"
    python3 /mnt/f/scripts/hic_scripts/diff_genes.py ${prefix}_diff.txt /mnt/f/organism_genome/Pfalciparum3D7/var.gff /mnt/f/organism_genome/Pfalciparum3D7/PlasmoDB-58_Pfalciparum3D7.chrom.sizes ${prefix}_diff.pdf
done
EOF
