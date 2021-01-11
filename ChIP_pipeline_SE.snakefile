import uuid
import os
from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
GS = GSRemoteProvider()

# Copy input bams from Google Bucket
rule get_input_files: 
	input: 
		GS.remote(expand("gs://landerlab-20210106-thouis-waszac-bams/{files}.bam", files=config["files"]))
		# "gs://landerlab-20210106-thouis-waszac-bams/{files}.bam"
		# "data/Waszak_2015/Waszak_2015_gs_names.txt"
	output: 
		temp("processed/{dataset}/aligned/{sample}.bam")
	resources:
		mem = 8000
	threads: 4
	shell:
		"""
		cp {input} {output}
		"""

#Sort BAM files by coordinates
rule sort_bams_by_position:
	input:
		"processed/{dataset}/aligned/{sample}.bam"
	output:
		protected("processed/{dataset}/sorted_bam/{sample}.sortedByCoords.bam")
	params:
		local_tmp = "/tmp/" + uuid.uuid4().hex + "/"
	resources:
		mem = 8000
	threads: 4
	shell:
		"""
		mkdir {params.local_tmp}
		cp {input} {params.local_tmp}/{wildcards.sample}.bam
		samtools sort -T {params.local_tmp}{wildcards.sample} -O bam -@ {threads} -o {params.local_tmp}{wildcards.sample}.sortedByCoords.bam {params.local_tmp}{wildcards.sample}.bam
		cp {params.local_tmp}{wildcards.sample}.sortedByCoords.bam {output}
		rm -r {params.local_tmp}
		"""
# samtools sort -T {params.local_tmp}/{wildcards.sample} -O bam -@ {threads} -o {params.local_tmp}/{wildcards.sample}.sortedByCoords.bam {params.local_tmp}/{wildcards.sample}.bam
# cp {params.local_tmp}/{wildcards.sample}.sortedByCoords.bam {output}

#Index sorted bams
rule index_bams:
	input:
		"processed/{dataset}/sorted_bam/{sample}.sortedByCoords.bam"
	output:
		"processed/{dataset}/sorted_bam/{sample}.sortedByCoords.bam.bai"
	resources:
		mem = 50
	threads: 1
	shell:
		"""
		samtools index {input}
		"""

#Keep only properly paired reads from nuclear chromosomes (remove MT)
rule filter_properly_paired:
	input:
		bam = "processed/{dataset}/sorted_bam/{sample}.sortedByCoords.bam",
		index = "processed/{dataset}/sorted_bam/{sample}.sortedByCoords.bam.bai"
	output:
		# temp("processed/{dataset}/filtered/{sample}.filtered.bam")
		"processed/{dataset}/filtered/{sample}.filtered.bam"
	resources:
		mem = 100
	threads: 1
	params:
		chr_list = "chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"
	shell:
		"""
		samtools view -h -b {input.bam} {params.chr_list} > {output}
		"""

#Remove bwa entry from the BAM file header (conflicts with MarkDuplicates)
rule remove_bwa_header:
	input:
		"processed/{dataset}/filtered/{sample}.filtered.bam"
	output:
		new_header = temp("processed/{dataset}/filtered/{sample}.new_header.txt"),
		bam = temp("processed/{dataset}/filtered/{sample}.reheadered.bam")
	resources:
		mem = 100
	threads: 1
	shell:
		"""
		# Remove bwa header
		samtools view -H {input} | grep -v 'ID:bwa' > {output.new_header}
		samtools reheader {output.new_header} {input} > {output.bam}
		"""

# Remove GI reference file header
rule remove_reference_header:
	input:
		"processed/{dataset}/filtered/{sample}.reheadered.bam"
	output:
		new_header = temp("processed/{dataset}/filtered/{sample}.new_header.txt"),
		bam = temp("processed/{dataset}/filtered/{sample}.rereheadered.bam")
	resources:
		mem = 100
	threads: 1
	shell:
		"""
		# Remove reference header
		samtools view -H {input} | grep -v 'SN:gi'
		samtools view -H {input} | grep -v 'SN:gi' > {output.new_header}
		samtools reheader {output.new_header} {input} > {output.bam}
		"""

# Reorder with Picard ReorderSam
rule reorder_bam:
	input:
		"processed/{dataset}/filtered/{sample}.rereheadered.bam"
	output:
		bam = "processed/{dataset}/filtered/{sample}.reordered.bam"
	resources:
		mem = 12000
	threads: 4
	shell:
		"""
		# samtools view -H {input}
		java -jar {config[picard_path]} ReorderSam INPUT={input} OUTPUT={output.bam} SEQUENCE_DICTIONARY={config[dictionary]} VALIDATION_STRINGENCY=LENIENT
		# samtools view -H {output}
		"""

#Remove duplicates using Picard MarkDuplicates
rule remove_duplicates:
	input:
		# "processed/{dataset}/filtered/{sample}.filtered.bam"
		# "processed/{dataset}/filtered/{sample}.reheadered.bam"
		"processed/{dataset}/filtered/{sample}.reordered.bam"
	output:
		bam = "processed/{dataset}/filtered/{sample}.no_duplicates.bam",
		metrics = "processed/{dataset}/metrics/{sample}.MarkDuplicates.txt"
		# bam = protected("processed/{dataset}/filtered/{sample}.no_duplicates.bam"),
		# metrics = protected("processed/{dataset}/metrics/{sample}.MarkDuplicates.txt")
	resources:
		mem = 12000
	threads: 4
	shell:
		"""
		# samtools view -H {input}
		java -jar {config[picard_path]} MarkDuplicates I={input} O={output.bam} REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT METRICS_FILE={output.metrics}
		"""

# Produce peaks
rule peaks:
    input:
        "processed/{dataset}/filtered/{sample}.no_duplicates.bam",
    output:
        narrowPeaks = protected("processed/{dataset}/peaks/{sample}.no_duplicates.narrowPeak"),
        summits = protected("processed/{dataset}/peaks/{sample}.no_duplicates.summits.bed"),
        saf_annot = protected("processed/{dataset}/peaks/{sample}.no_duplicates.saf"),
    resources:
        mem = 10000,
    threads:
        1,
    params:
        tmpOut = "/tmp/" + uuid.uuid4().hex + "/",
        tmpIn = "/tmp/" + uuid.uuid4().hex,
    	xls = "{sample}_peaks.xls",
    	tempPeaks = "{sample}_peaks.narrowPeak",
    	tempSummits = "{sample}_summits.bed",
    	saf_file = "processed/{dataset}/peaks/{sample}.no_duplicates.saf"
    shell:
    	"""
        mkdir {params.tmpOut} &&
        cp {input} {params.tmpIn} &&
        macs2 callpeak -t {params.tmpIn} -q 0.01 -n {wildcards.sample} --outdir {params.tmpOut} -f BAM --nomodel --extsize 147 &&
        # macs2 callpeak -t {params.tmpIn} -q 0.01 -n {wildcards.sample} --outdir {params.tmpOut} -f BAM &&
    	rm {params.tmpOut}{params.xls} &&
    	cp {params.tmpOut}{params.tempPeaks} {output.narrowPeaks} &&
        cp {params.tmpOut}{params.tempSummits} {output.summits} &&
        rm {params.tmpOut}{params.tempPeaks} &&
        rm {params.tmpOut}{params.tempSummits} &&
        rm -r {params.tmpOut} &&
        rm {params.tmpIn}
		# Convert .narrowPeak to saf file for featureCounts
		awk 'OFS="\t" {{print $1"."$2"."$3, $1, $2, $3, "."}}' {output.narrowPeaks} > {params.saf_file}
        """

rule feat_cnt:
    input:
        bam = "processed/{dataset}/filtered/{sample}.no_duplicates.bam",
        peak_annot = rules.peaks.output.saf_annot
    output:
        counts = "processed/{dataset}/counts/{sample}.no_duplicates.counts.txt",
        summary = "processed/{dataset}/counts/{sample}.no_duplicates.counts.txt.summary",
    threads:
        6,
    resources:
        mem = 20000,
    shell:
        """
		featureCounts -p -C -D 5000 -d 50 -F SAF -a {input.peak_annot} -o {output.counts} {input.bam}
        """
# featureCounts -p -C -D 5000 -d 50 -a {config[PU1_peaks_gtf]} -o {output.counts} {input.bam}
# featureCounts -p -C -D 5000 -d 50 -F GTF -t "exon" -a {input.peak_annot} -o {output.counts} {input.bam}

#Make sure that all final output files get created
rule make_all:
	input:
		#expand("processed/{{dataset}}/aligned/{sample}.bam", sample=config["samples"]),
		expand("processed/{{dataset}}/sorted_bam/{sample}.sortedByCoords.bam.bai", sample=config["samples"]),
		#expand("processed/{{dataset}}/filtered/{sample}.no_duplicates.bam", sample=config["samples"]),
		expand("processed/{{dataset}}/peaks/{sample}.no_duplicates.narrowPeak", sample=config["samples"]),
		expand("processed/{{dataset}}/counts/{sample}.no_duplicates.counts.txt", sample=config["samples"]),
	output:
		"processed/{dataset}/out.txt",
	resources:
		mem = 100,
	threads: 1,
	shell:
		"echo 'Done' > {output}"

# If input is fastq
#Align reads to the reference genome using bowtie2 aligner 
# rule align_with_bt_aln:
# 	input:
# 		"processed/{dataset}/fastq/{sample}.fastq.gz"
# 	output:
# 		temp("processed/{dataset}/aligned/{sample}.bam")
# 	params:
# 		index_dir = config["bt_index_dir"],
# 		index_name = config["bt_index_name"],
# 		basename = config["bt_basename"],
# 		rg="@RG\tID:{sample}\tSM:{sample}",
# 		tmp_fq = "/tmp/" + uuid.uuid4().hex + ".fastq.gz",
# 		tmp_sam = "/tmp/" + uuid.uuid4().hex + ".sam",
# 		tmp_bam = "/tmp/" + uuid.uuid4().hex + ".bam",
# 		tmp_index_dir = "/tmp/" + uuid.uuid4().hex + "/"
# 	resources:
# 		mem = 12000
# 	threads: 8
# 	shell:
# 		"""
# 		# head -10000 processed/PU1/fastq/NA06985_PU1.fastq.gz > processed/PU1/fastq/NA06985_PU1_small.fastq.gz
# 		cp {input} {params.tmp_fq}
# 		# cp processed/PU1/fastq/NA06985_PU1_small.fastq.gz {params.tmp_fq}		
# 		mkdir {params.tmp_index_dir}
# 		cp "{params.index_dir}"/* {params.tmp_index_dir}
# 		# ls {params.tmp_index_dir}
# 		bowtie2 -p {threads} -q --local -x {params.tmp_index_dir}{params.basename} -U {params.tmp_fq} -S {params.tmp_sam}
# 		samtools view -S -b {params.tmp_sam} > {params.tmp_bam}
# 		cp {params.tmp_bam} {output}
# 		rm {params.tmp_bam}
# 		rm {params.tmp_sam}
# 		rm {params.tmp_fq}
# 		rm -r {params.tmp_index_dir}
# 		"""
# bwa aln -t {threads} {params.tmp_index_dir}{params.index_name} {params.tmp_fq} > {params.tmp_sai}
# bwa samse -r '{params.rg}' {params.tmp_index_dir}{params.index_name} {params.tmp_sai} {params.tmp_fq} | samtools view -b - > {params.tmp_bam}

#Remove bowtie2 entry from the BAM file header (conflicts with MarkDuplicates)
# rule remove_bt_header:
# 	input:
# 		"processed/{dataset}/filtered/{sample}.filtered.bam"
# 	output:
# 		new_header = temp("processed/{dataset}/filtered/{sample}.new_header.txt"),
# 		bam = temp("processed/{dataset}/filtered/{sample}.reheadered.bam")
# 	resources:
# 		mem = 100
# 	threads: 1
# 	shell:
# 		"""
# 		# samtools view -H {input} | grep -v 'ID:bowtie2' > {output.new_header}
# 		# samtools view -H {input} | grep -v '^@' > {output.new_header}
# 		samtools reheader {output.new_header} {input} > {output.bam}
# 		"""
