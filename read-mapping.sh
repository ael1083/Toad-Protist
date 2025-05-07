#!/bin/bash

#Index reference genome
bwa index contigs.fasta

#Map reads to construct SAM file
bwa mem -t 24 contigs.fasta 16S_sequences.fasta > raw_mapped.sam
less -S raw_mapped.sam

#Construct a coverage table
samtools view -@ 24 -Sb  raw_mapped.sam  | samtools sort -@ 24 -o sorted_mapped.bam
samtools flagstat sorted_mapped.bam

#Index BAM file
samtools index sorted_mapped.bam

#Calculate per base coverage
bedtools genomecov -ibam sorted_mapped.bam > coverage.out
gen_input_table.py  --isbedfiles contigs.fasta coverage.out >  coverage_table.tsv
