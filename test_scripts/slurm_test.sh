#!/bin/bash
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 3-00:00:00
#SBATCH --mem-per-cpu=64G
#SBATCH -p intel

#SBATCH --mail-user=tlenz001@ucr.edu
#SBATCH --mail-type=end
#SBATCH --job-name=HiCpro_s1_todd_HiC
#SBATCH --export=ALL

mapfile -t avp_files < <(find /bigdata/lerochlab/shared/todd_HiC/toxo_project -name "*.allValidPairs")
for avp in "${avp_files[@]}"; do
    echo "removing sex chromosomes from ${avp}"
    python3 /mnt/f/scripts/test_scripts/test.py ${avp}
done

mapfile -t new_avp_files < <(find /bigdata/lerochlab/shared/todd_HiC/toxo_project -name "*_new.allValidPairs")
for new_avp in "${new_avp_files[@]}"; do
    python3 /mnt/f/scripts/hic_scripts/HiC-Pro/bin/utils/hicpro2juicebox.sh -i ${new_avp} -g /mnt/f/organism_genome/Hsapien_autosomes/GRCh38.chrom.sizes -j /mnt/f/scripts/hic_scripts/juicer/CPU/common/juicer_tools.jar -r /mnt/f/organism_genome/Hsapien_autosomes/GRCh38_MboI.bed
done
