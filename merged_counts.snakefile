import uuid
import os
from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
GS = GSRemoteProvider()

rule feat_cnt:
    input:
        bam = GS.remote("gs://landerlab-20210106-thouis-waszac-bams/PU1/{sample}.no_duplicates.bam"),
        peak_annot = "results/PU1/consensus.saf"
    output:
        counts = "processed/{dataset}/merged_counts/{sample}.merged_counts.txt",
        summary = "processed/{dataset}/merged_counts/{sample}.merged_counts.txt.summary",
    threads:
        6,
    resources:
        mem = 20000,
    shell:
        """
        featureCounts -p -C -D 5000 -d 50 -F SAF -a {input.peak_annot} -o {output.counts} {input.bam}
        """

#Make sure that all final output files get created
rule make_all:
	input:
		expand("processed/{{dataset}}/merged_counts/{sample}.merged_counts.txt", sample=config["samples"]),
	output:
		"processed/{dataset}/out.txt",
	resources:
		mem = 100,
	threads: 1,
	shell:
		"echo 'Done' > {output}"
