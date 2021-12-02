extract_raw_and_feature_fast.py: python script from Nanom6A version 18/03/2021.

predict_sites.py: python script from Nanom6A version 18/03/2021.

predict_sites_shuffled_WT_equal_KO.py: python script from Nanom6A version 18/03/2021 slightly modified to analyze mod and unmod filed shuffled from a WT dataset equal to the KO.

predict_sites_shuffled_WT_larger_KO.py: python script from Nanom6A version 18/03/2021 slightly modified to analyze mod and unmod filed shuffled from a WT dataset larger than the KO.

WT: folder containing Nanom6A results for the WT sample.

KD: folder containing Nanom6A results for the KO sample.

WT_shuffled: folder containing Nanom6A results for the WT sample with shuffled probabilities.

KD_shuffled: folder containing Nanom6A results for KO sample with shuffled probabilities.

WT_equal_KO: folder containing Nanom6A results for the WT sample sub-sampled to have a number of analyzed 5-mers equal to the KO.

WT_smaller_KO: folder containing Nanom6A results for the WT sample sub-sampled to have a number of analyzed 5-mers smaller than the KO.

WT_equal_shuffled: folder containing Nanom6A results for the WT sample sub-sampled with shuffled probabilities.

KO_equal_shuffled: folder containing Nanom6A results for the KO sample sub-sampled with shuffled probabilities.

GRCh37_latest_genomic.bed6.R: R script to create the bed6 file from RefSeq gtf file.

GRCh37_latest_genomic.fasta: genome fasta file from RefSeq. Not uploaded on GitHub due to its size, almost 1Gb after compression, it can be downloaded from: https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_genomic.fna.gz.

GRCh37_latest_genomic.fasta.fai: samtools index file for GRCh37_latest_genomic.fasta.

GRCh37_latest_genomic.dict: picard index file for GRCh37_latest_genomic.fasta.

GRCh37_latest_genomic.gtf.tar.gz: genome fasta file from RefSeq.

GRCh37_latest_genomic.bed6: bed file from RefSeq gtf.

refseq_GRCh37_latest_rna.fa.tar.gz: transcriptome fasta file from RefSeq.

refseq_GRCh37_latest_rna.fa.fai: samtools index file for refseq_GRCh37_latest_rna.fa.

refseq_GRCh37_latest_rna.dict: picard index file for refseq_GRCh37_latest_rna.fa.

probabilitiesShuffling.R: R script to generate mod and unmod files with shuffled probabilities.

wildTypeSubsampling.R: R script to sub-sample the WT sample.

dataAnalysis.R: R script to generate the comparative analysis shown in Figure 1.

comparativePlot.pdf: raw plot shown in Figure 1, generated through dataAnalysis.R.
