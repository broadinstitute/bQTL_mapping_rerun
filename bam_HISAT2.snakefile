from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
GS = GSRemoteProvider()

# Align reads to the reference genome using HISAT2 aligner 
rule align_with_hisat2_aln:
	input:
		GS.remote("gs://landerlab-20210106-thouis-waszac-bams/fastq/{sample}.fastq.gz"),
	output:
		bam = "processed/PU1/aligned_HISAT2/{sample}.bam",
		metrics = "processed/PU1/metrics_hisat/dup_metrics_{sample}.txt"
	resources:
		mem = 12000
	threads: 10
	shell:
		"""
		hisat2 --sensitive --no-discordant --no-spliced-alignment -p {threads} --remove-chrname \
			-x /mnt/DATA/grch37_snp/genome_snp \
			-U {input} | \
		bamsort inputformat=sam markduplicates=1 rmdup=1 fixmates=1 \
				inputthreads=8 outputthreads=8 M={output.metrics} \
				index=1 O={output.bam}
		"""

#Make sure that all final output files get created
rule make_all:
	input:
		expand("processed/PU1/aligned_HISAT2/{sample}.bam", sample=config["samples"]),
	output:
		"processed/PU1/out.txt",
	resources:
		mem = 100,
	threads: 1,
	shell:
		"echo 'Done' > {output}"

