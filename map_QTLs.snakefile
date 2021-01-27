from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
GS = GSRemoteProvider()

rule map_qtls:
	input:
		expand("processed/{{study}}/qtltools/output/{annot_type}/{condition}.permuted.txt.gz", annot_type = config["annot_type"], condition = config["conditions"]),
		expand("processed/{{study}}/qtltools/output/{annot_type}/sorted/{condition}.nominal.sorted.txt.gz", annot_type = config["annot_type"], condition = config["conditions"]),
		expand("processed/{{study}}/qtltools/output/{annot_type}/sorted/{condition}.nominal.sorted.txt.gz.tbi", annot_type = config["annot_type"], condition = config["conditions"]),
	output:
		"processed/{study}/out.txt"
	resources:
		mem = 100
	threads: 1
	shell:
		"echo 'Done!' > {output}"

#Compress and index input bed file
rule compress_bed:
	input:
		bed = "processed/{study}/qtltools/input/{annot_type}/{condition}.norm_prop.txt"
	output:
		bed = protected("processed/{study}/qtltools/input/{annot_type}/{condition}.norm_prop.txt.gz"),
		bed_index = protected("processed/{study}/qtltools/input/{annot_type}/{condition}.norm_prop.txt.gz.tbi")
	threads: 1
	resources:
		mem = 100
	shell:
		"""
		bgzip {input.bed} && tabix -p bed {output.bed}
		"""

#Run QTLtools in permutation mode
rule permutation_run:
	input:
		bed = "processed/{study}/qtltools/input/{annot_type}/{condition}.norm_prop.txt.gz",
		bed_index = "processed/{study}/qtltools/input/{annot_type}/{condition}.norm_prop.txt.gz.tbi",
		covariates = "processed/{study}/qtltools/input/{annot_type}/{condition}.covariates_prop.txt",
		# vcf = GS.remote("gs://landerlab-vcf/1000_genomes_vcfs/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"),
		# vcf_i = GS.remote("gs://landerlab-vcf/1000_genomes_vcfs/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi")
		vcf = "landerlab-vcf/1000_genomes_vcfs/Waszak_47_samples.chr1.vcf.gz"
	output:
		temp("processed/{study}/qtltools/output/{annot_type}/batches/{condition}.permutation.batch.{batch}.{n_batches}.txt")
	params:
		chunk = "{batch} {n_batches}"
	threads: 1
	resources:
		mem = 5000
	shell:
		"""
		# tabix -p vcf {input.vcf}
		QTLtools cis --vcf {input.vcf} --bed {input.bed} --cov {input.covariates} --chunk {params.chunk} --out {output} --window {config[cis_window]} --permute 10000 || touch {output}
		"""


#Merge all batches from QTLtools
rule merge_permutation_batches:
	input:
		expand("processed/{{study}}/qtltools/output/{{annot_type}}/batches/{{condition}}.permutation.batch.{batch}.{n_batches}.txt", 
			batch=[i for i in range(1, config["n_batches"] + 1)],
			n_batches = config["n_batches"])
	output:
		"processed/{study}/qtltools/output/{annot_type}/{condition}.permuted.txt.gz"
	resources:
		mem = 100
	threads: 1
	shell:
		"""
		cat {input} | bgzip > {output}
		"""


#Run QTLtools in nominal mode
rule nominal_run:
	input:
		bed = "processed/{study}/qtltools/input/{annot_type}/{condition}.norm_prop.txt.gz",
		bed_index = "processed/{study}/qtltools/input/{annot_type}/{condition}.norm_prop.txt.gz.tbi",
		covariates = "processed/{study}/qtltools/input/{annot_type}/{condition}.covariates_prop.txt",
		# vcf = GS.remote("gs://landerlab-vcf/1000_genomes_vcfs/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"),
		# vcf_i = GS.remote("gs://landerlab-vcf/1000_genomes_vcfs/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi")
		vcf = "landerlab-vcf/1000_genomes_vcfs/Waszak_47_samples.chr1.vcf.gz"
	output:
		temp("processed/{study}/qtltools/output/{annot_type}/nominal_batches/{condition}.nominal.batch.{batch}.{n_batches}.txt")
	params:
		chunk = "{batch} {n_batches}"
	threads: 1
	resources:
		mem = 5000
	shell:
		"""
		# tabix -p vcf landerlab-vcf/1000_genomes_vcfs/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
		QTLtools cis --vcf {input.vcf} --bed {input.bed} --cov {input.covariates} --chunk {params.chunk} --out {output} --window {config[nominal_cis_window]} --nominal 1 || touch {output}
		"""

#Merge all batches from QTLtools
rule merge_nominal_batches:
	input:
		expand("processed/{{study}}/qtltools/output/{{annot_type}}/nominal_batches/{{condition}}.nominal.batch.{batch}.{n_batches}.txt", 
			batch=[i for i in range(1, config["n_batches"] + 1)],
			n_batches = config["n_batches"])
	output:
		temp("processed/{study}/qtltools/output/{annot_type}/{condition}.nominal.txt.gz")
	resources:
		mem = 100
	threads: 1
	shell:
		"""
		cat {input} | bgzip > {output}
		"""

#Sort QTLtools output files
rule sort_qtltools_output:
	input:
		"processed/{study}/qtltools/output/{annot_type}/{condition}.nominal.txt.gz"
	output:
		protected("processed/{study}/qtltools/output/{annot_type}/sorted/{condition}.nominal.sorted.txt.gz")
	resources:
		mem = 1000
	threads: 2
	shell:
		"""
		zcat {input} | awk -v OFS='\\t' '{{$1=$1; print $0}}' | sort -k9,9 -k10,10n -k11,11n | bgzip > {output}
		"""

#Tabix-index QTLtools output files
rule index_qtltools_output:
	input:
		"processed/{study}/qtltools/output/{annot_type}/sorted/{condition}.nominal.sorted.txt.gz"
	output:
		"processed/{study}/qtltools/output/{annot_type}/sorted/{condition}.nominal.sorted.txt.gz.tbi"
	resources:
		mem = 1000
	threads: 1
	shell:
		"""
		tabix -s9 -b10 -e11 -f {input}
		"""

