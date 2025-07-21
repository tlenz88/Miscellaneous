#!/bin/bash
: '
# Set the sample basenames
conditions=("HFF" "ME49RFP_6hpi" "ME49RFP_24hpi")

# Set the resolution
resolutions=(1000000 500000 250000 100000 50000 10000)

for resolution in "${resolutions[@]}"; do 

    mapfile -t hic_files < <(find */juicer_files/ -name "*.hic")
    samples=()
    for file in "${hic_files[@]}"; do
        base_name=$(basename "${file%.hic}")
        samples+=("${base_name}")
    done

    echo Convert HiC to Cool format
    for sample in "${samples[@]}"; do
        hicConvertFormat -m "${sample}/juicer_files/${sample}.hic" -o "${sample}/hicexplorer_files/${resolution}/${sample}.cool" --inputFormat hic --outputFormat cool -r ${resolution}
    done

    mapfile -t cooler_files < <(find */hicexplorer_files/ -name "*_${resolution}.cool")
    sample_base=()
    h5_files=()
    norm_files=()
    for file in "${cooler_files[@]}"; do
        base_name="${file%.cool}"
        sample_base+=("${base_name}")
        h5_files+=("${base_name}.h5")
        norm_files+=("${base_name}_normalized.h5")
    done

    echo Convert Cool to H5 format
    for sample in "${sample_base[@]}"; do
        hicConvertFormat -m "${sample}.cool" -o "${sample}.h5" --inputFormat cool --outputFormat h5 -r "${resolution}"
    done
    
    echo Sum replicate matrices
    for condition in "${conditions[@]}"; do
        hicSumMatrices -m "${condition}_A_HiC/hicexplorer_files/${resolution}/${condition}_A_HiC_${resolution}.h5" \
        "${condition}_B_HiC/hicexplorer_files/${resolution}/${condition}_B_HiC_${resolution}.h5" \
        "${condition}_C_HiC/hicexplorer_files/${resolution}/${condition}_C_HiC_${resolution}.h5" \
        --outFileName "${condition}_HiC/hicexplorer_files/${resolution}/${condition}_HiC_${resolution}.h5"
    done

    merged_samples=("${conditions[0]}_HiC/hicexplorer_files/${resolution}/${conditions[0]}_HiC_${resolution}.h5" \
        "${conditions[1]}_HiC/hicexplorer_files/${resolution}/${conditions[1]}_HiC_${resolution}.h5" \
        "${conditions[2]}_HiC/hicexplorer_files/${resolution}/${conditions[2]}_HiC_${resolution}.h5")

    merged_normalized=()
    merged_sample_base=()
    for file in "${merged_samples[@]}"; do
        base_name="${file%.h5}"
        merged_sample_base+=("${base_name}")
        merged_normalized+=("${base_name}_normalized.h5")
    done

    echo Perform matrix normalization
    hicNormalize -m "${h5_files[@]}" -o "${norm_files[@]}" -n smallest
    hicNormalize -m "${merged_samples[@]}" -o "${merged_normalized[@]}" -n smallest

    for sample in "${sample_base[@]}"; do
        hicCorrectMatrix diagnostic_plot -m "${sample}_normalized.h5" -o "${sample}_normalized.png"
    done
    for sample in "${merged_sample_base[@]}"; do
        hicCorrectMatrix diagnostic_plot -m "${sample}_normalized.h5" -o "${sample}_normalized.png"
    done
    echo Perform matrix correction
    for sample in "${sample_base[@]}"; do
        hicCorrectMatrix correct -m "${sample}_normalized.h5" -o "${sample}_corrected.h5" --correctionMethod ICE -n 1000 --filterThreshold -4 3
    done
    for sample in "${merged_sample_base[@]}"; do
        hicCorrectMatrix correct -m "${sample}_normalized.h5" -o "${sample}_corrected.h5" --correctionMethod ICE -n 1000 --filterThreshold -4 3
    done

    echo Find TADs
    for sample in "${sample_base[@]}"; do
        hicFindTADs -m "${sample}_corrected.h5" --outPrefix "${sample}" --correctForMultipleTesting 'fdr' --threshold 0.05 -p 18
    done
    for sample in "${merged_sample_base[@]}"; do
        hicFindTADs -m "${sample}_corrected.h5" --outPrefix "${sample}" --correctForMultipleTesting 'fdr' --threshold 0.05 -p 18
    done

    echo Perform PCA
    for sample in "${sample_base[@]}"; do
        hicPCA --matrix "${sample}_corrected.h5" --format bigwig -o "${sample}_pca1.bw" "${sample}_pca2.bw" --whichEigenvectors 1 2 -p 18
    done
    for sample in "${merged_sample_base[@]}"; do
        hicPCA --matrix "${sample}_corrected.h5" --format bigwig -o "${sample}_pca1.bw" "${sample}_pca2.bw" --whichEigenvectors 1 2 -p 18
    done

    echo Generate A/B compartments
    for sample in "${sample_base[@]}"; do
        hicPCA --matrix "${sample}_corrected.h5" --format bedgraph -o "${sample}_compartments.bedgraph" --whichEigenvectors 1 -p 18
    done
    for sample in "${merged_sample_base[@]}"; do
        hicPCA --matrix "${sample}_corrected.h5" --format bedgraph -o "${sample}_compartments.bedgraph" --whichEigenvectors 1 -p 18
    done

    echo Detect loops
    for sample in "${sample_base[@]}"; do
        hicDetectLoops -m "${sample}_corrected.h5" -o "${sample}_loops.bedgraph" -p 18
    done
    for sample in "${merged_sample_base[@]}"; do
        hicDetectLoops -m "${sample}_corrected.h5" -o "${sample}_loops.bedgraph" -p 18
    done

    echo Transform to observed/expected matrix
    for sample in "${sample_base[@]}"; do
        hicTransform -m "${sample}_corrected.h5" --method obs_exp -o "${sample}_obs_exp.h5"
    done
    for sample in "${merged_sample_base[@]}"; do
        hicTransform -m "${sample}_corrected.h5" --method obs_exp -o "${sample}_obs_exp.h5"
    done

    echo Generate contact probability plots
    for sample in "${sample_base[@]}"; do
        hicPlotDistVsCounts -m "${sample}_corrected.h5" -o "${sample}_distance_counts.png"
    done
    for sample in "${merged_sample_base[@]}"; do
        hicPlotDistVsCounts -m "${sample}_corrected.h5" -o "${sample}_distance_counts.png"
    done

    echo Compare two samples
    hicCompareMatrices -m "${merged_sample_base[1]}".h5 "${merged_sample_base[0]}".h5 -o "${conditions[1]}_vs_${conditions[0]}".h5
    hicCompareMatrices -m "${merged_sample_base[2]}".h5 "${merged_sample_base[1]}".h5 -o "${conditions[2]}_vs_${conditions[1]}".h5

done
: '
while IFS=$'\t' read -r chr len; do
    for i in WT_subsample/hicexplorer_files/5000; do
        hicPlotTADs --tracks ${i}/${i}_${resolution}_tracks.ini --outFileName ${i}/${i}_${resolution}_${chr}.pdf --region ${chr}:0-${len} --fontSize 8 --width 40
    done
done < /mnt/f/organism_genome/CparvumBGF/CpBGF_genome_v16.chrom.sizes

