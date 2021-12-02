#!/bin/bash

## CPUs
cpus=number_of_CPUs

## Setup references
genomeBed=GRCh37_latest_genomic.bed6

genomeFasta=GRCh37_latest_genomic.fasta

## Setup Nanom6A scripts
predict_sites=predict_sites.py

## Setup Nanom6A models
models=/path/to/nanom6A/models

## Setup required folders
nanom6A_step1=/path/to/nanom6A/step1/results/WT_equal_KO/

nanom6A_step2=/path/to/nanom6A/step2/results/WT_equal_KO/

## Nanom6A step 2
cd $nanom6A_step2
conda activate nanom6A_step2
python $predict_sites --model $models --cpu $cpus -i $nanom6A_step1/result -o result_final -r $genomeBed -g $genomeFasta --support 20
conda deactivate