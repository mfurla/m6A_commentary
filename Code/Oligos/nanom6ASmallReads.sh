#!/bin/bash

## CPUs
cpus=number_of_CPUs

## Setup required folders
fast5sTmp=/folder/for/fast5s/involved/in/the/analysis

fastqsTmp=/folder/for/fastqs/involved/in/the/analysis

filenames=/path/to/folder/for/filenames/

nanom6A_step1=/path/to/nanom6A/step1/results/

nanom6A_step2=/path/to/nanom6A/step2/results/

tombo=/path/to/tombo/folder/

## Setup references
transcriptomeFasta=oligo_reference.fa
transcriptomeBed=oligo_reference.bed

## Setup Nanom6A scripts
extract_raw_and_feature_fast=extract_raw_and_feature_fast.py
predict_sites_smallReads=predict_sites_smallReads.py

## Setup Nanom6A models
models=/path/to/nanom6A/models

## Setup data
fast5files=/path/to/all/passed/fast5s/*.fast5

## Install Tombo
conda activate tombo
cd $tombo
pip install -e .
conda deactivate

## Selection of long reads and multi to single read fast5s conversion
for filename in $fast5files; do

	echo $filename

	# Temporary variables required to work on fastq and fast5 files
	tmp=$(echo "${filename##*/}")
	tmp=$(echo "${tmp%%.*}")

	# Retrieve the reads longer than 100 bps for the analysis from the fastq file
	cd $fastqsTmp

	# From fast5 name to fastq
	cp ${filename//fast5/fastq} .

	# Assign to each read a length and select those longer than 100 bps
	grep "runid" $tmp".fastq" | cut -c2-37 > readsIds.txt 
	grep -A1 "runid" $tmp".fastq" | grep -v "\-\-" | grep -v "@" >  reads.txt
	awk '{ print length }' reads.txt > readsLength.txt
	
	paste readsIds.txt readsLength.txt > readsAndLength.txt
	
	awk '{if ($2<100) {print $1".fast5"} }' readsAndLength.txt > shortReadsIds.txt

	cp readsAndLength.txt readsAndLength_$tmp.txt

	# Extract the single read fast5 files and select the long ones
	cd $fast5sTmp

	echo "FILE COPY: $filename"
	cp $filename .

	# Multi to single reads fast5s
	conda activate ont-fast5-api
	multi_to_single_fast5 --input_path ${tmp}.fast5 --save_path .
	rm -f FAL*
	mv 0/* .
	rm -r 0
	conda deactivate

	# Remove the files annotated in the shortReadsIds.txt file
	xargs rm < $fastqsTmp/shortReadsIds.txt

done

## Tombo resquiggle	
cd $fast5sTmp

conda activate tombo
tombo resquiggle . $transcriptomeFasta --basecall-group Basecall_1D_000 --overwrite --processes $cpus --fit-global-scale --include-event-stdev --failed-reads-filename "failedReads_$tmp.txt"
conda deactivate
 
## Remove reads badly processed by Tombo
grep "Poor raw to expected signal matching" $"failedReads_$tmp.txt" > poorSignalReads.txt
sed -i 's/Poor raw to expected signal matching (revert with `tombo filter clear_filters`)//g' poorSignalReads.txt
sed -i 's/BaseCalled_template:::.\//\n/g' poorSignalReads.txt
sed -i 's/,//g' poorSignalReads.txt
		
xargs rm < $fast5sTmp/poorSignalReads.txt

## Remove reads not mapped
grep "Alignment not produced" $"failedReads_$tmp.txt" > notAlignedReads.txt
sed -i 's/Alignment not produced//g' notAlignedReads.txt
sed -i 's/BaseCalled_template:::.\//\n/g' notAlignedReads.txt
sed -i 's/,//g' notAlignedReads.txt
		
xargs rm < $fast5sTmp/notAlignedReads.txt

mv $"failedReads_$tmp.txt" $filenames

## Nanom6A step 1
cd $nanom6A_step1

find $fast5sTmp -maxdepth 1 -name "*.fast5" > files.txt

conda activate nanom6A_step1
python $extract_raw_and_feature_fast --cpu=$cpus  --fl=/path/to/files.txt -o "result_$tmp" --clip=0 # The m6A site is close to the end of the reads
conda deactivate

cat *.tsv > result.feature.tsv
cat *.fa > result.feature.fa

## Nanom6A step 2
cd $nanom6A_step2

conda activate nanom6A_step2
python $predict_sites_smallReads --model $models --cpu $cpus  -i $nanom6A_step1/result -o $nanom6A_step2 -r $transcriptomeBed -g $transcriptomeFasta --proba 0.5
conda deactivate