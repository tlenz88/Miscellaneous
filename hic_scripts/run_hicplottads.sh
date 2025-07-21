while IFS=$'\t' read -r chr len; do 
	for i in 500000 250000 100000 50000; do 
		for j in HFF ME49RFP_6hpi ME49RFP_24hpi; do
			hicPlotTADs \
				--tracks /mnt/f/toxo_project/HiC/Hsapien_output/output_files/${j}_HiC/hicexplorer_files_two_reps/${i}/${j}_${i}.ini \
				--outFileName /mnt/f/toxo_project/HiC/Hsapien_output/output_files/${j}_HiC/hicexplorer_files_two_reps/${i}/${j}_${i}_${chr}_tracks.pdf \
				--region ${chr}:0-${len}
		done
	done
done < /mnt/f/organism_genome/Hsapien_autosomes/GRCh38.chrom.sizes
