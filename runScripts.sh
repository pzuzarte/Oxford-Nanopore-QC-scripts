#!/bin/bash

source activate py311

cd ~/Desktop/Oxford-Nanopore-QC-scripts

python ont_qc.py \
    --file ~/Desktop/FukDam/sequencing_summary_PBK65018_0067cf0f_8aa78e9d.txt \
    --runName FukDam \
    --outdir ~/Desktop/FukDam/QC_report/ \
    --maxLength 1000000 \
    --prop 100000 \
    --logLength \
    --video


python ont_qc.py \
    --file ~/Desktop/HelKat/sequencing_summary_PBK65121_66e5f4e7_a476615f.txt \
    --runName HelKat \
    --outdir ~/Desktop/HelKat/QC_report/ \
    --maxLength 1000000 \
    --prop 100000 \
    --logLength \
    --video
