: <<"COMMENT"
for i in Astro HFF; do
	for j in D3 D5 D7; do
		bamCompare \
		-b1 "${i}_${j}_K4_rep1/${i}_${j}_K4_rep1_dedup.bam" \
		-b2 "${i}_${j}_K4_rep2/${i}_${j}_K4_rep2_dedup.bam" \
		-o "${i}_${j}_K4.bw" \
		-of bigwig \
		--operation mean \
		--scaleFactorsMethod None \
		--normalizeUsing CPM \
		--ignoreDuplicates \
		-bs 1 -p 12
	done
done

for i in Astro HFF; do 
	for j in D3 D5 D7; do
		computeMatrix reference-point \
		-R "ToxoDB-62_TgondiiME49.gtf" \
		-S "${i}_${j}_K4.bw" \
		-o "${i}_${j}_K4_matrix.gz" \
		-a 1000 -b 1000 \
		-bs 10 -p 12
	done
done

for i in Astro HFF; do
	for j in D3 D5 D7; do
		plotHeatmap -m "${i}_${j}_K4_matrix.gz" \
		-o "${i}_${j}_K4_heatmap.pdf" \
		--colorMap Reds \
		--whatToShow 'heatmap and colorbar' \
		--outFileSortedRegions "${i}_${j}.bed" \
		--outFileNameMatrix "${i}_${j}_matrix.bed.gz"
	done
done

for i in Astro HFF; do 
	for j in D3 D5 D7; do
		computeMatrix scale-regions \
		-R "ToxoDB-62_TgondiiME49.gtf" \
		-S "${i}_${j}_K4.bw" \
		-o "${i}_${j}_K4_matrix.gz" \
		-bs 10 -p 12 -b 1000
	done
done
COMMENT
for i in Astro HFF; do
	for j in D3 D5 D7; do
		plotHeatmap -m "all_genes_coding_region/${i}_${j}_K4_matrix.gz" \
		-o "all_genes_coding_region/${i}_${j}_K4_heatmap.pdf" \
		--colorMap Reds \
		--whatToShow 'heatmap and colorbar' \
		--outFileSortedRegions "all_genes_coding_region/${i}_${j}.bed" \
		--outFileNameMatrix "all_genes_coding_region/${i}_${j}_matrix.bed.gz" \
		-max 15.3671
	done
done

for i in Astro HFF; do
	for j in D3 D5 D7; do
		plotHeatmap -m "all_genes_TSS_and_coding_region/${i}_${j}_K4_matrix.gz" \
		-o "all_genes_TSS_and_coding_region/${i}_${j}_K4_heatmap.pdf" \
		--colorMap Reds \
		--whatToShow 'heatmap and colorbar' \
		--outFileSortedRegions "all_genes_TSS_and_coding_region/${i}_${j}.bed" \
		--outFileNameMatrix "all_genes_TSS_and_coding_region/${i}_${j}_matrix.bed.gz" \
		-max 16.4042
	done
done

for i in Astro HFF; do
	for j in D3 D5 D7; do
		plotHeatmap -m "all_genes_TSS/${i}_${j}_K4_matrix.gz" \
		-o "all_genes_TSS/${i}_${j}_K4_heatmap.pdf" \
		--colorMap Reds \
		--whatToShow 'heatmap and colorbar' \
		--outFileSortedRegions "all_genes_TSS/${i}_${j}.bed" \
		--outFileNameMatrix "all_genes_TSS/${i}_${j}_matrix.bed.gz" \
		-max 17.2810
	done
done