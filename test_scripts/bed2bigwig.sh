#!/bin/bash

bamCoverage -b ${1} \
-o ${1%.*}.bw \
-bs 1 \
--normalizeUsing CPM \
-p 8 \
--effectiveGenomeSize 2875001522
