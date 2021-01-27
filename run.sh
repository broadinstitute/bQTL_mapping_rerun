#Process PU.1 data
#Convert bams to fastq files
# snakemake processed/PU1/out.txt -s bam_to_fastq_SE.snakefile --configfile configs/config_PU1.yaml --cores 5

# Login for gcloud for google bucket access
# gcloud auth login

# #Run the processing pipeline
# snakemake -n processed/PU1/out.txt -s ChIP_pipeline_SE.snakefile -j 10 --configfile configs/config_PU1.yaml --cores 5

# Extract genotypes
# snakemake -s munge/PU1/extract_genotypes.snakefile munge/PU1/out.txt --configfile configs/config_PU1.yaml --cores 5

#Map QTLs
snakemake -s map_QTLs.snakefile processed/PU1/out.txt -j 10 --configfile configs/config_PU1.yaml  --cores 5