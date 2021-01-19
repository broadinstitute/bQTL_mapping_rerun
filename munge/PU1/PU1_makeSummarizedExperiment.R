library("dplyr")
library("cqn")
# library("devtools")
library("SummarizedExperiment")
library("rtracklayer")
library("seqUtils")

# Set working directory
setwd('/mnt/DATA/Blood_ATAC')

# Import Peaks and GC-content 
peaks = readr::read_delim("results/PU1/GCcontent.txt", delim = "\t", col_types = "ciiddiiiiiii") 
colnames(peaks) = c("chr","start","end","pct_at","percentage_gc_content","num_A","num_C","numG","num_T","num_N", "num_oth","length")
regions = peaks[,c(1,2,3,12)] 
regions = regions %>% as.data.frame() %>% tibble::as_tibble() %>% dplyr::mutate(strand = 1)
gc_content = peaks[,c(1,2,3,5)]

# Import phenotype IDs
phenotype_IDs = readr::read_delim("results/PU1/consensus.saf", delim = "\t", col_names = FALSE, col_types = "cciic") 
colnames(phenotype_IDs) = c("phenotype_id","chr","start","end","strand")
phenotype_IDs = phenotype_IDs[,c(1,2,3,4)]

# Merge & label row names
peak_data = dplyr::left_join(regions, gc_content, by = c("chr","start","end")) %>% as.data.frame() 
peak_data = dplyr::left_join(phenotype_IDs, peak_data, by = c("chr","start","end")) %>% as.data.frame() 
rownames(peak_data) = phenotype_IDs$phenotype_id
peak_data = dplyr::rename(peak_data, gene_id = phenotype_id)

#Import counts
counts = readr::read_delim("results/PU1/counts.merged.txt",  delim = "\t")
count_matrix = dplyr::select(counts, -Chr,-Start,-End,-Strand, -Geneid, -Length) %>% as.matrix()
rownames(count_matrix) = phenotype_IDs$phenotype_id
genotype_id = c("NA06985_PU1","NA06986_PU1","NA06994_PU1","NA07037_PU1","NA07048_PU1","NA07051_PU1","NA07056_PU1","NA07346_PU1","NA07357_PU1","NA10847_PU1","NA10851_PU1","NA11829_PU1","NA11830_PU1","NA11831_PU1","NA11832_PU1","NA11840_PU1","NA11881_PU1","NA11894_PU1","NA11918_PU1","NA11920_PU1","NA11931_PU1","NA11992_PU1","NA11993_PU1","NA11994_PU1","NA12005_PU1","NA12043_PU1","NA12154_PU1","NA12156_PU1","NA12234_PU1","NA12249_PU1","NA12275_PU1","NA12282_PU1","NA12286_PU1","NA12287_PU1","NA12383_PU1","NA12489_PU1","NA12750_PU1","NA12760_PU1","NA12761_PU1","NA12762_PU1","NA12763_PU1","NA12776_PU1","NA12812_PU1","NA12813_PU1","NA12814_PU1","NA12815_PU1","NA12873_PU1")
colnames(count_matrix) = genotype_id

# Normalize with CQN
cqn_matrix = calculateCQN(count_matrix, peak_data)
# head(cqn_matrix)

#Import the list of genotyped samples
genotyped = read.table("data/Waszak_2015/Waszak_2015_genotype_list.txt")

#Import 1000G metadata
kg_meta = readr::read_delim("metadata/1000_genomes_sample_metadata.txt",  delim = "\t")
kg_meta = dplyr::rename(kg_meta, genotype_id = Sample)
# head(kg_meta)

 #Import Waszak metadata
sample_meta = readr::read_tsv("data/Waszak_2015/Waszak_2015_clean_metadata.txt") %>% dplyr::filter(sample_id %in% colnames(count_matrix)) %>% dplyr::select(sample_id, genotype_id) %>% dplyr::left_join(kg_meta, by = "genotype_id") %>% dplyr::filter(genotype_id %in% genotyped$V1) %>% as.data.frame()
rownames(sample_meta) = sample_meta$sample_id
# head(sample_meta)

#Construct a SummarizedExperiment object
se = SummarizedExperiment::SummarizedExperiment(
  assays = list(counts = count_matrix[,rownames(sample_meta)], cqn = cqn_matrix[,rownames(sample_meta)]), 
  colData = sample_meta, 
  rowData = peak_data)
saveRDS(se, "results/SummarizedExperiment/PU1_SummarizedExperiment.rds")