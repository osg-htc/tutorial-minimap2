#!/bin/bash

# Use minimap2 to map the basecalled reads to the reference genome
# The command takes three arguments: reference genome, reads file, and output SAM file
# Example usage: ./run_minimap2.sh <ref_genome> <reads.fastq> <output_sam_file>

# Define arguments/variables
ref_genome="$1"
reads_file="$2"
output_sam_file="$3"

# Run minimap2 with ONT preset (-x map-ont)
minimap2 -ax map-ont "$ref_genome" "$reads_file" > "$output_sam_file"
echo "Minimap2 mapping complete. Output SAM file: $output_sam_file"

# Use samtools to sort our mapped reads BAM, required for downstream analysis
output_bam_file="${output_sam_file%.sam}_sorted.bam"
samtools sort "$output_sam_file" -o "$output_bam_file"
echo "Samtools sorting complete. Output sorted BAM file: $output_bam_file"
