# Order of files for data analysis 

### 1. merged_counts.sh
* Concatenate each sample's count files into one matrix of peaks x samples
* Environment: r_4.0.3

### 2. makeSummarizedExperiment.ipynb
* Merge count data and sample metadata info
* Environment: r_4.0.3

### 3. export_QTLtools.ipynb
* Preprocess and prepare matrix for QTLtools 
* Environment: r_4.0.3

### 4. extract_genotypes.snakefile
* Extract the 45 samples from VCF
* Environment: PU1_env

### 5. merge_QTLtools_PCA.sh
* PCA of phenotype (and genotype, optional)
* Environment: qtltools
* <code> bash merge_QTLtools_PCA.sh </code>

### 6. merge_QTLtools_PCA.ipynb
* Construct covariate matrix with phenotype and genotype PCAs
* Environment: r_4.0.3

### 7. map_qtls.snakefile
* Map bQTLs
* This is in the main directory
* Environment: PU1_env

### 8. QTL_check.ipynb
* Produce list of significant bQTLs
* Environment: r_4.0.3

### 9. QTL_to_bed.sh (optional)
* Convert list of significant bQTLs to .bed file (ex. to view on IGV)
* <code> bash QTL_to_bed.sh </code>

## Visualizations
### 10. Visualizations_PCA_Reads_Peaks.ipynb
* Visualizations of count matrix and PCAs (phenotype and genotype)
* Environment: PU1_env

### 11. Compare_PCA_Methods.ipynb
* Comparing genotype PCA using 45 samples, CEU subpopulation, or full 1KG population
* Environment: PU1_env

## Other
### 12. merged_consensus_peaks.sh
* Merges the 47 peak files into 1 files of consensus peaks 
* This is integrated ChIP_pipeline_SE.snakefile as rule consensus_peaks but this script produces a histogram relevant to deciding threshold for "consensus"
