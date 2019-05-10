#!/bin/bash

input_dir=$1
target_bed=$2
bedtools_genome_file=$3

mapping_qc() {
    
    # run samtools to get statistics on mapped reads
    echo '## Flagstat before filtering off-target reads' > "$1"".stats"
    samtools flagstat "$1" >> "$1"".stats"
    
    # run bedtools to determine on-target reads
    intersectBed -a "$1" -b $target_bed -u -wa > "$1"".ontarget.bam"
    
    # run samtools to get statistics on the on-target mapped reads
    echo '## Flagstat after filtering off-target reads' >> "$1"".stats"
    samtools flagstat "$1"".ontarget.bam" >> "$1"".stats"

    # clean up 
    rm "$1"".ontarget.bam"

    # run bedtools to generate a histogram coverage files
    coverageBed -hist -a $target_bed -b "$1" -sorted -g $bedtools_genome_file | grep ^all > "$1"".coverage.report"
    
    # run bedtools to generate coverage for target bed
    coverageBed -a $target_bed -b "$1" -sorted -g $bedtools_genome_file > "$1"".coverage.report.bed"

}

for bam in "$input_dir"/*/*.bam; do mapping_qc "$bam" & done


