#!/bin/bash

echo "SORTING GFF FILE BY CHROMOSOME AND COORDINATE"
sort -t $'\t' -k1,1 -k4,4n /mnt/d/organism_genome/Pfalciparum3D7/PlasmoDB-50_${1}.gff > /mnt/d/organism_genome/Pfalciparum3D7/gff_chr/PlasmoDB-50_${1}sorted.gff

echo "SPLITTING GFF FILE BY CHROMOSOME"
python3 /mnt/d/scripts/split_gff_by_chr.py /mnt/d/organism_genome/Pfalciparum3D7/gff_chr/PlasmoDB-50_${1}sorted.gff /mnt/d/organism_genome/Pfalciparum3D7/gff_chr/

echo "SPLITTING MeDIPseq DATA BY CHROMOSOME"
for i in ./*; do echo $i; mkdir ${i}/bed_chr; python3 /mnt/d/scripts/split_chr.py ${i}/*_sorted_cpm.bed ${i}/bed_chr/${i}; done

echo "SPLITTING MeDIPseq DATA BY CHROMOSOME"
for i in {01..14}; do for j in *; do python /mnt/d/scripts/medip_genes.py ./${j}/bed_chr/${j}_Pf3D7_${i}_v3.txt /mnt/d/organism_genome/Pfalciparum3D7/gff_chr/Pf3D7_${i}_v3.gff ./${j}/bed_chr/${j}_Pf3D7_${i}_v3_genecov.txt; done; done

echo "MERGING FEATURE COUNTS FILES"d
for i in *; do cat ./${1}/bed_chr/*_genecov.txt > ./${1}/${1}_genecov.txt; done