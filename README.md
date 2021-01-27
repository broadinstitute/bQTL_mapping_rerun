# Reanalysis of ChIP-seq data for TF PU.1 in lymphoblastoid cell lines

This repostiory contains reanalysis of ChIP-seq data from the following paper:

* Waszak, S. M. et al. Population variation and genetic control of modular chromatin architecture in humans. *Cell.* **162**, 1039–1050 (2015).

The inputs are bams of the 47 samples for PU.1 and the outputs are the signficant binding QTLs (bQTLs) for chromosome 1. 

Fastq files can also be inputs (uncomment rules in the bottom of ChIP_pipeline_SE.snakefile to align with bowtie2).

Credits - Fork of https://github.com/kauralasoo/Blood_ATAC. I only looked at Waszak et al 2015 data.