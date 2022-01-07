# CHROMS = ["1","10","11","12","13","14","15","16","17","18","19","2","20","21","22","3","4","5","6","7","8","9","X","Y"]
CHROMS = ["3","4","5","6","7","8","9","10"]

rule make_all:
    input:
        # expand("/gpfs/hpchome/a72094/rocket/projects/chromatin-QTLs/bQTL_mapping_rerun/results/PU1/genotypes/Waszak_47_samples.chr{chromosome}.vcf.gz", chromosome = CHROMS)
        expand("/home/jupyter/DATA/1000_genomes_vcfs/Waszak_47_samples.chr{chromosome}.vcf.gz", chromosome = CHROMS)
    output:
        "munge/{study}/out.txt"
    threads: 1
    resources:
        mem = 1000
    shell:
        "echo 'Done! > {output}'" 
        
rule extract_genotypes:
    input:
        # vcf = "/gpfs/hpchome/a72094/rocket/datasets/1000G/GRCh38/ALL.chr{chromosome}_GRCh38.genotypes.20170504.vcf.gz",
        # vcf = GS.remote("gs://landerlab-vcf/1000_genomes_vcfs/ALL.{chromosome}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"),
        vcf = "/home/jupyter/DATA/1000_genomes_vcfs/ALL.chr{chromosome}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",
        samples = "/home/jupyter/DATA/bQTL_mapping_rerun/data/Waszak_2015/Waszak_2015_genotype_list.txt"
    output:
        vcf = "/home/jupyter/DATA/1000_genomes_vcfs/Waszak_47_samples.chr{chromosome}.vcf.gz"
    threads: 1
    resources:
        mem = 3000
    shell:
        """
        bcftools view -S {input.samples}  --force-samples {input.vcf} | bcftools filter -i 'MAF[0] >= 0.05' -O z - > {output.vcf}
        """

#Merge all samples
#bcftools concat Waszak_47_samples.chr1.vcf.gz Waszak_47_samples.chr2.vcf.gz Waszak_47_samples.chr10.vcf.gz Waszak_47_samples.chr11.vcf.gz Waszak_47_samples.chr12.vcf.gz Waszak_47_samples.chr13.vcf.gz Waszak_47_samples.chr14.vcf.gz Waszak_47_samples.chr15.vcf.gz Waszak_47_samples.chr16.vcf.gz Waszak_47_samples.chr17.vcf.gz Waszak_47_samples.chr18.vcf.gz Waszak_47_samples.chr19.vcf.gz Waszak_47_samples.chr20.vcf.gz Waszak_47_samples.chr21.vcf.gz Waszak_47_samples.chr22.vcf.gz Waszak_47_samples.chr3.vcf.gz Waszak_47_samples.chr4.vcf.gz Waszak_47_samples.chr5.vcf.gz Waszak_47_samples.chr6.vcf.gz Waszak_47_samples.chr7.vcf.gz Waszak_47_samples.chr8.vcf.gz Waszak_47_samples.chr9.vcf.gz -O z -o Waszak_47_samples.merged.vcf.gz


#Merge chr1-10 samples
#bcftools concat Waszak_47_samples.chr1.vcf.gz Waszak_47_samples.chr10.vcf.gz Waszak_47_samples.chr2.vcf.gz Waszak_47_samples.chr3.vcf.gz Waszak_47_samples.chr4.vcf.gz Waszak_47_samples.chr5.vcf.gz Waszak_47_samples.chr6.vcf.gz Waszak_47_samples.chr7.vcf.gz Waszak_47_samples.chr8.vcf.gz Waszak_47_samples.chr9.vcf.gz -O z -o Waszak_47_samples.merged_chr1_10.vcf.gz