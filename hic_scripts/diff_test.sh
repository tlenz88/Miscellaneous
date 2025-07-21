
gene_file="/mnt/f/toxo_project/HiC/Hsapien_output/output_files/genes_of_interest.txt"
resolution=10000

while read -r gene chrom start end; do
    plot_start=$(( ((start + resolution / 2) / resolution * resolution) - 1000000 ))
    if [[ $plot_start -lt 0 ]]; then
        plot_start=0
    fi
    plot_end=$(( ((start + resolution / 2) / resolution * resolution) + 1000000 ))
    echo "${chrom}:${plot_start}-${plot_end}"
    hicPlotTADs --tracks "/mnt/f/toxo_project/HiC/Hsapien_output/output_files/diff_${resolution}_${chrom}.ini" \
    --region "${chrom}:${plot_start}-${plot_end}" \
    --outFileName "/mnt/f/toxo_project/HiC/Hsapien_output/output_files/${gene}_region.pdf"
done < "$gene_file"
