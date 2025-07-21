#!/bin/bash

# List the resolutions, sample names, and replicates
resolutions=(10000)
samples=("HFF" "ME49RFP_6hpi" "ME49RFP_24hpi")

for res in "${resolutions[@]}"; do
    for sample in "${samples[@]}"; do
        # Set basename for input/output files
        bname="${sample}_HiC/hicexplorer_files_two_reps/${res}/${sample}_${res}"
        : '
        # Plot histogram of bin coverage to identify matrix correction threshold
        hicCorrectMatrix diagnostic_plot \
        -m "${bname}_normalized.h5" \
        -o "${bname}_normalized.png" \
        --perchr \
        &> "${bname}_normalized_mad.txt" 
        : '
        # Extract threshold for matrix correction
        for chrom in 4 8 9 15 16 17 18 20; do 
            case "$chrom" in
                4|8|15|20)
                    if [[ ${sample} == "HFF" ]]; then
                        continue
                    fi
                ;;
            esac
            case "$chrom" in
                4|8|9|16|18)
                    if [[ ${sample} == "ME49RFP_6hpi" ]]; then
                        continue
                    fi
                ;;
            esac
            case "$chrom" in
                16|17|20)
                    if [[ ${sample} == "ME49RFP_24hpi" ]]; then
                        continue
                    fi
                ;;
            esac
            echo
            echo
            echo
            echo $sample $chrom
            # Extract the mad threshold for the current chromosome
            madscore=$(grep "mad threshold" "${bname}_normalized_mad2.txt" | \
            grep "chr${chrom}:" | \
            sed -E "s/.*chr${chrom}: mad threshold ([^ ]+).*/\1/")

            # Check if madscore is empty or contains invalid data
            if [[ -z "$madscore" || ! "$madscore" =~ ^-?[0-9]*\.?[0-9]+$ ]]; then
                echo "Error: Invalid MAD threshold for chr${chrom}."
                continue
            fi

            madscore=-3
            # Calculate the upper threshold as -3 times the madscore
            upper=$(echo "-3 * $madscore" | bc)

            # Run hicCorrectMatrix using the threshold
            hicCorrectMatrix correct \
            -m "${bname}_normalized.h5" \
            -o "${bname}_chr${chrom}_corrected.h5" \
            --correctionMethod ICE \
            --filterThreshold $madscore $upper \
            --chromosomes chr${chrom} \
            --skipDiagonal

            # Identify TADs
            hicFindTADs -m "${bname}_chr${chrom}_corrected.h5" \
            --outPrefix "${bname}_chr${chrom}" \
            --correctForMultipleTesting 'fdr' \
            --thresholdComparisons 0.05 \
            -p 12
            : '
            # Plot read count over distance
            hicPlotDistVsCounts -m "${bname}_chr${chrom}_corrected.h5" \
            -o "${bname}_chr${chrom}_dvc.pdf" \
            -s
            : '
            if [[ ${res} -ge 50000 ]]; then
                if [[ "${sample}" == "HFF" ]]; then
                    extraTrackRep="F"
                else
                    extraTrackRep="D"
                fi

                # Identify A/B compartments
                hicPCA --matrix "${bname}_chr${chrom}_corrected.h5" \
                -o "${bname}_chr${chrom}_PCA1.bw" "${bname}_chr${chrom}_PCA2.bw" \
                --whichEigenvectors 1 2 \
                -oem "${bname}_chr${chrom}_oem.h5" \
                -pm "${bname}_chr${chrom}_pm.h5" \
                --extraTrack "/mnt/f/toxo_project/CUTandTag/Hsapien_output/output_files/${sample}_H3K27ac_${extraTrackRep}_CT/${sample}_H3K27ac_${extraTrackRep}_CT_filtered.bw" \
                --histonMarkType "active" \
                -p 12
                
                hicTransform --matrix "${bname}_chr${chrom}_corrected.h5" \
                --outFileName "${bname}_chr${chrom}_pearson.h5" \
                --method pearson \
                --perChromosome

                # Identify loops
                hicDetectLoops -m "${bname}_chr${chrom}_corrected.h5" \
                -o "${bname}_chr${chrom}_loops.bedgraph"
            fi
        done
    done

    for (( i=0; i<${#samples[@]}; i++ )); do
        for (( j=0; j<${#samples[@]}; j++ )); do
            if [ "$i" -eq "$j" ]; then
                continue
            fi

            diff_out="${samples[i]}_HiC/hicexplorer_files_two_reps/${res}"

            for chrom in 4 8 9 15 16 17 18 20; do
                # Find differential interactions between samples
                hicCompareMatrices -m "${diff_out}/${samples[i]}_${res}_chr${chrom}_corrected.h5" \
                "${samples[j]}_HiC/hicexplorer_files_two_reps/${res}/${samples[j]}_${res}_chr${chrom}_corrected.h5" \
                -o "${diff_out}/${samples[i]}_vs_${samples[j]}_${res}_chr${chrom}_normalized.h5" \
                --operation log2ratio

                # Find differential TADs between samples
                hicDifferentialTAD -tm "${diff_out}/${samples[i]}_${res}_chr${chrom}_corrected.h5" \
                -cm "${samples[j]}_HiC/hicexplorer_files_two_reps/${res}/${samples[j]}_${res}_chr${chrom}_corrected.h5" \
                -td "${diff_out}/${samples[i]}_${res}_chr${chrom}_domains.bed" \
                -o "${diff_out}/${samples[i]}_vs_${samples[j]}_${res}_chr${chrom}_diff_TADs_one" \
                -p 0.05 \
                -m all \
                -mr one \
                -t 12

                # Find differential TADs between samples
                hicDifferentialTAD -tm "${diff_out}/${samples[i]}_${res}_chr${chrom}_corrected.h5" \
                -cm "${samples[j]}_HiC/hicexplorer_files_two_reps/${res}/${samples[j]}_${res}_chr${chrom}_corrected.h5" \
                -td "${diff_out}/${samples[i]}_${res}_chr${chrom}_domains.bed" \
                -o "${diff_out}/${samples[i]}_vs_${samples[j]}_${res}_chr${chrom}_diff_TADs_all" \
                -p 0.05 \
                -m all \
                -mr all \
                -t 12

                rename 's/.diff_tad/.bed/' "${diff_out}/${samples[i]}_vs_${samples[j]}_${res}_chr${chrom}_diff_TADs_one_accepted.diff_tad"
                rename 's/.diff_tad/.bed/' "${diff_out}/${samples[i]}_vs_${samples[j]}_${res}_chr${chrom}_diff_TADs_all_accepted.diff_tad"
                rename 's/.diff_tad/.bed/' "${diff_out}/${samples[i]}_vs_${samples[j]}_${res}_chr${chrom}_diff_TADs_one_rejected.diff_tad"
                rename 's/.diff_tad/.bed/' "${diff_out}/${samples[i]}_vs_${samples[j]}_${res}_chr${chrom}_diff_TADs_all_rejected.diff_tad"
            done
        done
    done
done
