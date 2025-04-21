#!/bin/bash

#Run PROKKA to do genome annotations
nohup prokka --centre X --compliant contigs.fasta --outdir prokka_output --cpus 24 --mincontiglen 200 &
ls prokka_output

#Get counts for each gene annotation

