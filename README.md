# Toad-Protist

## Background
Sequencing Platform
- Paired-End
- 250 bp sequencing reads

These files include:
- 1 total sample

Goals:
- Assemble genome of protist
- Identify protist
- Assess genome

## Methods
### tool
- tools used and what they do

## Results

## Code
```bash
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

#Run QUAST to provide basic statistics
quast.py contigs.fasta -o quast_results

#Run BUSCO to assess completeness of genome assembly
busco -i contigs.fasta -m genome -o busco-results -l bacteria

#Run PROKKA to do genome annotations
nohup prokka --centre X --compliant contigs.fasta --outdir prokka_output --cpus 24 --mincontiglen 200 &
ls prokka_output

#Get counts for each gene annotation
grep -o "product=.*" prokka_output/PROKKA_*.gff | sed 's/product=//g' | sort | uniq -c | sort -nr > protein_abundances.txt

#Extract 16S sequences from FFN file
filter_fasta_by_taxonomy_and_length.py --keys "16S ribosomal RNA" --out 16S_sequences.fasta prokka_output/PROKKA_*.ffn

#Make BLAST database out of contig assembly
makeblastdb -in contigs.fasta -dbtype nucl -out contigs_db

#Run BLAST
blastn -query 16S_sequences.fasta -db contigs_db -out 16S_vs_contigs_6.tsv -o
utfmt 6
less 16S_vs_contigs_6.tsv

#Create output file for blobtools
/usr/local/bin/blast-ncbi-nt.sh contigs.fasta
less contigs.fasta.vs.nt.cul5.1e5.megablast.out
```
