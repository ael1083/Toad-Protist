#!/bin/blast

#Make BLAST database out of contig assembly
makeblastdb -in contigs.fasta -dbtype nucl -out contigs_db

#Run BLAST
blastn -query 16S_sequences.fasta -db contigs_db -out 16S_vs_contigs_6.tsv -o
utfmt 6
less 16S_vs_contigs_6.tsv

#Create output file for blobtools
/usr/local/bin/blast-ncbi-nt.sh contigs.fasta
less contigs.fasta.vs.nt.cul5.1e5.megablast.out
