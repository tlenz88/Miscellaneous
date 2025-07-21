metadata="$1"

Rscript /mnt/f/scripts/dcHiC/dchicf.r --file "$metadata" --pcatype cis --dirovwt T
Rscript /mnt/f/scripts/dcHiC/dchicf.r --file "$metadata" --pcatype select --dirovwt T --genome hg38
Rscript /mnt/f/scripts/dcHiC/dchicf.r --file "$metadata" --pcatype analyze --dirovwt T
Rscript /mnt/f/scripts/dcHiC/dchicf.r --file "$metadata" --pcatype viz --genome hg38

Rscript /mnt/f/scripts/dcHiC/dchicf.r --file "$metadata" --pcatype subcomp --dirovwt T
Rscript /mnt/f/scripts/dcHiC/dchicf.r --file "$metadata" --pcatype fithic --dirovwt T --fithicpath "/mnt/f/scripts/hic_scripts/fithic/fithic/fithic.py" --pythonpath "~/miniconda3/lib/python3.12"
Rscript /mnt/f/scripts/dcHiC/dchicf.r --file "$metadata" --pcatype dloop --dirovwt T
Rscript /mnt/f/scripts/dcHiC/dchicf.r --file "$metadata" --pcatype viz --genome hg38
Rscript /mnt/f/scripts/dcHiC/dchicf.r --file "$metadata" --pcatype enrich --genome hg38