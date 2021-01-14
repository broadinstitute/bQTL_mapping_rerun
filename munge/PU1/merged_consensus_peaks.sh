#!/bin/bash

# Parameters
INPUT_PATH="processed/PU1/peaks"
OUTPUT_PATH="results/PU1"
REF_GENOME="landerlab-bowtie-hg19-index-misc/genome.fa"
THRESHOLD=30

# Empty result directory
rm $OUTPUT_PATH/*

# Concatenate all narrowPeak peaks (sample.narrowPeak --> allpeaks.bed)
cat ${INPUT_PATH}/*.narrowPeak > ${OUTPUT_PATH}/allpeaks.bed

# Sort merged bed file (allpeaks.bed --> allpeaks_sorted.bed)
bedtools sort -i ${OUTPUT_PATH}/allpeaks.bed > ${OUTPUT_PATH}/allpeaks_sorted.bed

# Produce counts (genome coverage track)
bedtools coverage -a ${OUTPUT_PATH}/allpeaks_sorted.bed -b ${OUTPUT_PATH}/allpeaks_sorted.bed -sorted -counts > ${OUTPUT_PATH}/allpeaks_counts.bed

# Produce histogram of counts to decide filter threshold
# bedtools coverage -a ${OUTPUT_PATH}/allpeaks_sorted.bed -b ${OUTPUT_PATH}/allpeaks_sorted.bed -sorted -hist > ${OUTPUT_PATH}/allpeaks_counts_hist.txt
awk '{ print $NF }' ${OUTPUT_PATH}/allpeaks_counts.bed > ${OUTPUT_PATH}/allpeaks_counts.txt
gnuplot << \EOF
set style histogram
set xrange [0:48]
set style fill solid 1.0 noborder
width=1
hist(x,width)=width/2.0 + width*floor(x/width)
set term png
set output 'results/PU1/histogram.png'
plot 'results/PU1/allpeaks_counts.txt' using (hist(column(1),width)):(1.0) smooth freq w boxes 
EOF

# Filter: extract just those lines where the last column is at least 23 > goodpeaks.bed
awk -v threshold=$THRESHOLD '$NF>=threshold' ${OUTPUT_PATH}/allpeaks_counts.bed > ${OUTPUT_PATH}/goodpeaks.bed

# Output consensus peaks
bedtools merge -i ${OUTPUT_PATH}/goodpeaks.bed > ${OUTPUT_PATH}/consensus.bed

# Output GC content (%GC content in column 5)
bedtools nuc -fi $REF_GENOME -bed ${OUTPUT_PATH}/consensus.bed > ${OUTPUT_PATH}/GCcontent.txt

# Done
echo "Done"
