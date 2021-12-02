#!/bin/bash

## CPUs
cpus=number_of_CPUs

## Setup required folders
FAST5_trash=/path/for/fast5s/not/used

nanom6A_step1=/path/to/nanom6A/step1/results/WT/

nanom6A_step2=/path/to/nanom6A/step2/results/WT/

## Setup references
transcriptomeFasta=refseq_GRCh37_latest_rna.fa

genomeBed=GRCh37_latest_genomic.bed6

genomeFasta=GRCh37_latest_genomic.fasta

## Setup Nanom6A scripts
extract_raw_and_feature_fast=extract_raw_and_feature_fast.py
predict_sites=predict_sites.py

## Setup Nanom6A models
models=/path/to/nanom6A/models

## Setup data
fast5files=/path/to/all/passed/fast5s/*.fast5

## Working directory
cd /path/to/results/folder/

## Resquiggle
conda activate tombo
tombo resquiggle $fast5files --ignore-read-locks --overwrite --basecall-group Basecall_1D_001 $transcriptomeFasta --processes $cpus --fit-global-scale --include-event-stdev --failed-reads-filename "failedReads.txt"
conda deactivate

## Problematic reads removal
# Bad tombo's performance
conda activate ont-fast5-api
grep -o '[^, ]\+' failedReads.txt | grep -o -P '(?<=BaseCalled_template:::).*' > problematicReads.txt
cat problematicReads.txt | parallel -j $cpus mv {} $FAST5_trash

# Alignment to reverse strand
find $fast5files -maxdepth 2 -name "*.fast5" | parallel -j $cpus h5ls -rvd {} | grep -A3 -e 'mapped_strand' -e 'fast5' >> tomboStrandsTmp.txt

grep -B 8 '"-"' tomboStrandsTmp.txt | grep -o -P '(?<=Opened ").*(?=" with sec2 driver.)' > tomboStrands.txt
cat tomboStrands.txt | parallel -j $cpus mv {} $FAST5_trash
conda deactivate

## Final FAST5s
find $fast5files -maxdepth 2 -name "*.fast5" > files.txt

## Nanom6A step 1
cd $nanom6A_step1
conda activate nanom6A_step1
python $extract_raw_and_feature_fast --cpu=$cpus --fl=/path/to/files.txt -o result --clip=10
conda deactivate

## Nanom6A step 2
cd $nanom6A_step2
conda activate nanom6A_step2
python $predict_sites --model $models --cpu $cpus -i $nanom6A_step1/result -o result_final -r $genomeBed -g $genomeFasta --support 20
conda deactivate