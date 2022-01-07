# Process PU.1 data

# Login for gcloud for google bucket access
# gcloud auth login

#Convert bams to fastq files
# snakemake processed/PU1/out.txt -s bam_to_fastq_SE.snakefile -j 10 --configfile configs/config_PU1.yaml --cores 5

# #Run the processing pipeline
snakemake processed/PU1/out.txt -s ChIP_pipeline_SE.snakefile -j 10 --configfile configs/config_PU1.yaml --cores 5

# Extract genotypes
# snakemake munge/PU1/out.txt -s munge/PU1/extract_genotypes.snakefile --configfile configs/config_PU1.yaml --cores 5

#Map QTLs
# snakemake -s map_QTLs.snakefile processed/PU1/out.txt -j 10 --configfile configs/config_PU1.yaml  --cores 5

# Align with HISAT2
# snakemake -s bam_HISAT2.snakefile processed/PU1/out.txt -j 10 --configfile configs/config_PU1.yaml  --cores 5