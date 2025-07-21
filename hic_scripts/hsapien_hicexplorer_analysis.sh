#!/bin/bash

# List the resolutions, sample names, and replicates
resolutions=(100000)
samples=("HFF" "ME49RFP_6hpi" "ME49RFP_24hpi")
#replicates=("A" "B" "C")

for res in "${resolutions[@]}"; do
    all_merged=() # Array for merged replicates
    all_merged_normalized=() # Array for normalized merged replicates
    all_merged_samples=() # Array for merged sample names
    all_merged_corrected=() # Array for corrected merged replicates
    all_merged_domains=() # Array for TADs
    : '
    all_replicates=() # Array for individual replicates
    all_replicates_normalized=() # Array for normalized replicates
    all_replicate_samples=() # Array for replicate sample names
    all_replicates_corrected=() # Array for corrected replicates
    : '
    for sample in "${samples[@]}"; do
        : '
        # Convert .hic file to .cool for given resolution
        hic_file="${sample}_HiC/juicer_files/${sample}_HiC_autosomes.hic"

        hicConvertFormat -m "${hic_file}" \
        -o "${sample}_HiC/hicexplorer_test/${res}/${sample}.cool" \
        --inputFormat hic \
        --outputFormat cool \
        -r "$res"
        : '
        bname="${sample}_HiC/hicexplorer_test/${res}/${sample}_${res}"
        : '
        # Convert .cool file to .h5 for given resolution
        hicConvertFormat -m "${bname}.cool" \
        -o "${bname}.h5" \
        --inputFormat cool \
        --outputFormat h5 \
        -r "$res"
        : '
        # Add file and sample names to arrays for use in normalization
        all_merged+=("${bname}.h5")
        all_merged_normalized+=("${bname}_normalized.h5")
        all_merged_samples+=("${sample}")
        : '
        for rep in "${replicates[@]}"; do
            rep_hic_file="${sample}_${rep}_HiC/juicer_files/${sample}_${rep}_HiC_autosomes.hic"

            hicConvertFormat -m "${rep_hic_file}" \
            -o "${sample}_${rep}_HiC/hicexplorer_test/${res}/${sample}_${rep}.cool" \
            --inputFormat hic \
            --outputFormat cool \
            -r "$res"

            rbname="${sample}_${rep}_HiC/hicexplorer_test/${res}/${sample}_${rep}_${res}"

            hicConvertFormat -m "${rbname}.cool" \
            -o "${rbname}.h5" \
            --inputFormat cool \
            --outputFormat h5 \
            -r "$res"

            all_replicates+=("${rbname}.h5")
            all_replicates_normalized+=("${rbname}_normalized.h5")
            all_replicate_samples+=("${sample}_${rep}")
        done
        : '
    done
    # Normalize all matrices to the lowest read count among samples
    hicNormalize -m "${all_merged[@]}" \
    -o "${all_merged_normalized[@]}" \
    -n smallest
    : '
    hicNormalize -m "${all_replicates[@]}" \
    -o "${all_replicates_normalized[@]}" \
    -n smallest

    # Plot sample correlation between all replicates at given resolution
    hicCorrelate -m "${all_replicates_normalized[@]}" \
    --method pearson \
    --log1p \
    --labels "${all_replicate_samples[@]}" \
    -oh "correlation_heatmap_${res}.pdf" \
    -os "correlation_scatterplot_${res}.pdf" \
    --plotNumbers
    : '
    for sample in "${samples[@]}"; do
        # Set basename for input/output files
        bname="${sample}_HiC/hicexplorer_test/${res}/${sample}_${res}"

        # Plot histogram of bin coverage to identify matrix correction threshold
        hicCorrectMatrix diagnostic_plot \
        -m "${bname}_normalized.h5" \
        -o "${bname}_normalized.png" \
        --perchr \
        &> "${bname}_normalized_mad.txt" 

        # Extract threshold for matrix correction
        madscore=$(grep "mad threshold" \
        "${bname}_normalized_mad.txt" | \
        sed 's/INFO:hicexplorer.hicCorrectMatrix:mad threshold //g')
        upper=$(echo -3*$madscore | bc)

        # Use threshold to perform iterative correction on matrices
        hicCorrectMatrix correct \
        -m "${bname}_normalized.h5" \
        -o "${bname}_corrected2.h5" \
        --correctionMethod ICE \
        --filterThreshold $madscore $upper \
        --perchr

        all_merged_corrected+=("${bname}_corrected.h5")

        # Identify TADs
        hicFindTADs -m "${bname}_corrected.h5" \
        --outPrefix "${bname}" \
        --correctForMultipleTesting 'fdr' \
        --thresholdComparisons 0.05 \
        -p 12

        all_merged_domains+=("${bname}_domains.bed")

        # Plot read count over distance
        hicPlotDistVsCounts -m "${bname}_corrected.h5" \
        -o "${bname}_dvc.png"

        if [[ ${res} -ge 50000 ]]; then
            if [[ "${sample}" == "HFF" ]]; then
                extraTrackRep="F"
            else
                extraTrackRep="D"
            fi

            # Identify A/B compartments
            hicPCA --matrix "${bname}_corrected.h5" \
            -o "${bname}_PCA1.bw" "${bname}_PCA2.bw" \
            --whichEigenvectors 1 2 \
            -oem "${bname}_oem.h5" \
            -pm "${bname}_pm.h5" \
            --extraTrack /mnt/f/toxo_project/CUTandTag/Hsapien_output/output_files/${sample}_H3K27ac_${extraTrackRep}_CT/${sample}_H3K27ac_${extraTrackRep}_CT_filtered.bw \
            --histonMarkType "active" \
            -p 12

            # Identify loops
            hicDetectLoops -m "${bname}_corrected.h5" \
            -o "${bname}_loops.bedgraph"
        fi
        : '
        for rep in "${replicates[@]}"; do
            rbname="${sample}_${rep}_HiC/hicexplorer_test/${res}/${sample}_${rep}_${res}"

            hicCorrectMatrix diagnostic_plot \
            -m "${rbname}_normalized.h5" \
            -o "${rbname}_normalized.png" \
            &> "${rbname}_normalized_mad.txt"

            replicate_madscore=$(grep "mad threshold" \
            "${rbname}_normalized_mad.txt" | \
            sed 's/INFO:hicexplorer.hicCorrectMatrix:mad threshold //g')
            replicate_upper=$(echo -3*$replicate_madscore | bc)

            hicCorrectMatrix correct \
            -m "${rbname}_normalized.h5" \
            -o "${rbname}_corrected.h5" \
            --correctionMethod ICE \
            --filterThreshold $replicate_madscore $replicate_upper

            all_replicates_corrected+=("${rbname}_corrected.h5")

            hicFindTADs -m "${rbname}_corrected.h5" \
            --outPrefix "${rbname}" \
            --correctForMultipleTesting 'fdr' \
            --thresholdComparisons 0.05 \
            -p 12

            hicPlotDistVsCounts -m "${rbname}_corrected.h5" \
            -o "${rbname}_dvc.png"

            if [[ ${res} -ge 50000 ]]; then
                if [[ "${sample}" == "HFF" ]]; then
                    extraTrackRep="F"
                else
                    extraTrackRep="D"
                fi

                hicPCA --matrix "${rbname}_corrected.h5" \
                -o "${rbname}_PCA1.bw" "${rbname}_PCA2.bw" \
                --whichEigenvectors 1 2 \
                -oem "${rbname}_oem.h5" \
                -pm "${rbname}_pm.h5" \
                --extraTrack /mnt/f/toxo_project/CUTandTag/Hsapien_output/output_files/${sample}_H3K27ac_${extraTrackRep}_CT/${sample}_H3K27ac_${extraTrackRep}_CT_filtered.bw \
                --histonMarkType "active" \
                -p 12

                # Identify loops
                hicDetectLoops -m "${rbname}_corrected.h5" \
                -o "${rbname}_loops.bedgraph"
            fi
        done
        : '
    done
    : '
    # Plot ratio of short vs long range interactions
    hicPlotSVL -m "${all_merged_corrected[@]}" \
    -pfn "merged_short_vs_long_range_interactions_${res}.pdf" \
    -o "merged_short_vs_long_range_interactions_${res}.txt" \
    -t 12

    hicPlotSVL -m "${all_replicates_corrected[@]}" \
    -pfn "replicate_short_vs_long_range_interactions_${res}.pdf" \
    -o "replicate_short_vs_long_range_interactions_${res}.txt" \
    -t 12

    for (( i=0; i<${#samples[@]}; i++ )); do
        for (( j=0; j<${#samples[@]}; j++ )); do
            if [ "$i" -eq "$j" ]; then
                continue
            fi

            diff_out="${all_merged_samples[i]}_HiC/hicexplorer_test/${res}"

            # Find differential interactions between samples
            hicCompareMatrices -m "${all_merged_corrected[i]}" \
            "${all_merged_corrected[j]}" \
            -o "${diff_out}/${all_merged_samples[i]}_vs_${all_merged_samples[j]}_${res}_normalized.h5" \
            --operation log2ratio

            # Find differential TADs between samples
            hicDifferentialTAD -tm "${all_merged_corrected[i]}" \
            -cm "${all_merged_corrected[j]}" \
            -td "${all_merged_domains[i]}" \
            -o "${diff_out}/${all_merged_samples[i]}_vs_${all_merged_samples[j]}_${res}_diff_TADs_one" \
            -p 0.05 \
            -m all \
            -mr one \
            -t 12

            # Find differential TADs between samples
            hicDifferentialTAD -tm "${all_merged_corrected[i]}" \
            -cm "${all_merged_corrected[j]}" \
            -td "${all_merged_domains[i]}" \
            -o "${diff_out}/${all_merged_samples[i]}_vs_${all_merged_samples[j]}_${res}_diff_TADs_all" \
            -p 0.05 \
            -m all \
            -mr all \
            -t 12

            rename 's/.diff_tad/.bed/' "${diff_out}/${samples[i]}_vs_${samples[j]}_${res}_diff_TADs_one_accepted.diff_tad"
            rename 's/.diff_tad/.bed/' "${diff_out}/${samples[i]}_vs_${samples[j]}_${res}_diff_TADs_all_accepted.diff_tad"
            rename 's/.diff_tad/.bed/' "${diff_out}/${samples[i]}_vs_${samples[j]}_${res}_diff_TADs_one_rejected.diff_tad"
            rename 's/.diff_tad/.bed/' "${diff_out}/${samples[i]}_vs_${samples[j]}_${res}_diff_TADs_all_rejected.diff_tad"
        done
    done
    : '
done
