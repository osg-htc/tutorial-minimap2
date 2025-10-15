#!/bin/bash

# Use minimap2 to map the basecalled reads to the reference genome
# The command takes two arguments: reference genome and reads file
# Example usage: ./run_minimap2.sh <ref_genome> <reads.fastq>

# Define arguments/variables
ref_genome=$1
reads_file=$2

# Run minimap2 with ONT preset (-x map-ont)
minimap2 -ax map-ont ${ref_genome} ${reads_file} > mapped_${reads_file}_reads_to_genome_sam
echo "Minimap2 mapping complete. Output SAM file: mapped_${reads_file}_reads_to_genome_sam"

# Use samtools to sort our mapped reads BAM, required for downstream analysis
output_bam_file="mapped_${reads_file}_reads_to_genome_sam_sorted.bam"
samtools sort "mapped_${reads_file}_reads_to_genome_sam" -o ${output_bam_file}
echo "Samtools sorting complete. Output sorted BAM file: ${output_bam_file}"