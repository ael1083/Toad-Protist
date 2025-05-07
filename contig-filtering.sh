#!/bin/bash

#Create lookup table
blobtools create -i contigs.fasta -b sorted_mapped.bam -t contigs.fasta.vs.nt.cul5.maxt10.1e5.megablast.out -o blob_out

#Create output table & plot
blobtools view -i blob_out.blobDB.json -r all -o blob_taxonomy
grep -v '##' blob_taxonomy.blob_out.blobDB.table.txt
blobtools plot -i blob_out.blobDB.json -r genus

#Filter by length (<500bp) & coverage (>0)
mkdir mdibl-t3-2018-WGS/filtered_assembly
cp ../blob_taxonomy.blob_out.blobDB.table.txt ./
grep -v '#' blob_taxonomy.blob_out.blobDB.table.txt | awk -F'\t' '$2 < 5000' | awk -F'\t' '$5 > 0' | awk -F'\t' '{print $1}' > list_of_contigs_to_keep_len500_cov20.txt
less -S list_of_contigs_to_keep_len5000_cov<0.txt

#Filter assembly based on contigs list
filter_contigs_by_list.py ~/toad-analysis/trimmed-reads/spades_assembly_default/contigs.fasta list_of_contigs_to_keep_len500_cov20.txt protist_filtered.fasta

#Find average coverage
grep -f list_of_contigs_to_keep_len500_cov20.txt blob_taxonomy.blob_out.blobDB.table.txt | awk '{w = w + $2; e = e + $5 * $2;} END {print e/w}'

#BLAST final contigs
wget "https://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec"
blastn -reward 1 -penalty -5 -gapopen 3 -gapextend 3 -dust yes -soft_masking true -evalue 700 -searchsp 1750000000000 -query protist_filtered.fasta -subject UniVec  -outfmt 6 -out genome_vs_univec.6
