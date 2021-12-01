resquiggle.py: slightly modified version of the same file from Tombo v1.5.1. Specifically, I set a seed for numpy.random, and I changed minimap2 parameters to facilitate the alignment of short reads.

nanom6ASmallReads.sh: bash script to run nanom6A analysis on the oligos datasets.

oligo_reference.fa: fasta file for the synthetic oligos.

oligo_reference.fa.fai: samtools index file for oligo_reference.fa.

oligo_reference.dict: picard index file for oligo_reference.fa.

oligo_reference.bed: fasta file for the synthetic oligos.

extract_raw_and_feature_fast.py: python script from Nanom6A version 11/06/2020.

predict_sites_smallReads.py: python script from Nanom6A version 11/06/2020 slightly modified to facilitate the alignment of short reads and manipulate the modification probability threshold (for exploratory studies).

control: folder containing Nanom6A results for the control oligo.

oligo1: folder containing Nanom6A results for the oligo with m6A.

oligosPerformance.R: script R to estimate Nanom6A classification performance.

rocCurve.pdf: analysis of Nanom6A classification performance, result of oligosPerformance.R.
