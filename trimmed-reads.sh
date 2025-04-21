#!/bin/bash

#Ran Fastqc on paired samples
mkdir fastqc_trimmed-reads
fastqc Fecal_S68_L001_R1_001.fastq.gz unpaired-Fecal_S68_L001.fastq.gz -o fastqc_trimmed-reads

#Run SPAdes to assemble the genomes
nohup spades.py -1 Fecal_S68_L001_R1_001.fastq.gz -2 Fecal_S68_L001_R2_001.fastq.gz -s unpaired-Fecal_S68_L001_R1_001.fastq.gz -s unpaired-Fecal_S68_L001_R2_001.fastq.gz -o spades_assembly_default -t 24 &

#View output data
grep ">" spades_assembly_default/contigs.fasta | head
grep -c '>' spades_assembly_default/contigs.fasta

#Clean
cd spades_spades_assembly_default/
mv contigs.fasta spades.log ../
rm -r *
mv ../contigs.fasta ../spades.log ./
ls

