import uuid
import os
from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
GS = GSRemoteProvider()

#Convert BAM files to fastq
rule bam_to_fastq:
	input:
		# "processed/{dataset}/input_bam/{sample}.bam"
		GS.remote("gs://landerlab-20210106-thouis-waszac-bams/{sample}.bam")
	output:
		"processed/{dataset}/fastq/{sample}.fastq.gz"
	resources:
		mem = 1000
	threads: 1
	shell:
		"samtools fastq -F 2816 -0 {output} -c 6 {input}"

#Make sure that all final output files get created
rule make_all:
	input:
		expand("processed/{{dataset}}/fastq/{sample}.fastq.gz", sample=config["samples"])
	output:
		"processed/{dataset}/out.txt"
	resources:
		mem = 100
	threads: 1
	shell:
		"echo 'Done' > {output}"

#module load samtools-1.6


