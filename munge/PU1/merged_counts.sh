#!/bin/bash

# Merge 47 count files into one matrix of peaks (2662) x samples (47)

# Parameters
INPUT_PATH="processed/PU1/merged_counts"
TEMP_FILE="results/PU1/TEMP.txt"
OUTPUT_TEMP_FILE="results/PU1/TEMP_OUTPUT.txt"
OUTPUT_FILE="results/PU1/counts.merged.txt"

rm $OUTPUT_FILE

first=true
for f in ${INPUT_PATH}/*.txt; do
    if [ $first == true ]
    then
        # Copy the geneid, chr, start, end, strand, length, and counts columns from first file
        tail -n +2 $f > $OUTPUT_TEMP_FILE 
        first=false
    else 
        # Copy only the last column (counts) from all other files
        awk '{print $NF}' $f | tail -n +2 > $TEMP_FILE
        paste $OUTPUT_TEMP_FILE $TEMP_FILE > $OUTPUT_FILE
        cp -f $OUTPUT_FILE $OUTPUT_TEMP_FILE
    fi
done

rm $TEMP_FILE
rm $OUTPUT_TEMP_FILE