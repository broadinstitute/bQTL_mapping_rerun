#!/bin/bash

# In order to view on IGV, change QTL list to a bed file

# File paths
INPUT_FILE="/home/jupyter/DATA/bQTL_mapping_rerun/processed/PU1/qtltools/output/cqn/results.genes.significant.txt"
TEMP_FILE="/home/jupyter/DATA/bQTL_mapping_rerun/processed/PU1/qtltools/output/cqn/TEMP.txt"
OUTPUT_PATH="/home/jupyter/DATA/bQTL_mapping_rerun/processed/PU1/qtltools/output/cqn/results.genes.significant.bed"

rm ${TEMP_FILE}

# Change delimiter to tab
sed 's/ /\t/g' ${INPUT_FILE} > ${TEMP_FILE}

# Remove first column
cut -f 2- ${TEMP_FILE} > ${OUTPUT_PATH}

rm ${TEMP_FILE}
