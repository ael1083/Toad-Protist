# Toad-Protist Project

Alexandria Lyden, Lawrence Gordon, W. Kelley Thomas

<details> <summary><H2> Background </H2></summary>

The data for this analysis was provided by the Hubbard Center for Genome Studies at the University of New Hampshire. It consists of 16s data in paired-end 250 bp reads that were amplified by Illumina HiSeq 2500. The files were made up of 1 sample taken from the Anaxyrus americanus, known as the American Toad. With this data, the goal was to assemble the genome of the protist, identify the protist's class, and assess the genome. The protist, Amphibiothecum penneri, was identified from [this](https://pubmed.ncbi.nlm.nih.gov/16456158/) paper.

</details></details>

<details> <summary><H2> Methods </H2></summary>

### Trimmomatic
This tool is used to trim low quality bases determined from a fastqc output. It will also remove adaptors automatically. This tool was used in conjunction with a script that made the program easier to use. The input is the raw forward and reverse reads and the outputs are trimmed fastq files.

### SPADes
This tool assembles bacterial genomes. The inputs are the trimmed fastq files and the output is a fasta file containing the genome assembly.

### QUAST
This is a tool that examines how well the genome assembly has been constructed and gives the total genome size. It was used to keep track of the length of each contig and provided statistics about them. The input is the genome assembly fasta file and the output consist of multiple tables.

### BUSCO
This tool assesses the completeness of the genome assembly using the OrthoDB set of single-copy orthologous that are found in at least 90% of all the organisms. The input was the genome assembly of contigs and the bacteria database and the output is a directory containing a summary of the results, a table with coordinates for where each orthologous gene is located in the assembly, and a directory with the nucleotide and amino acid sequences of all the identified sequences.

### PROKKA
This tool is used to annotate the genome assembly. The inputs are the contig fasta file, the output is gene annotations in GFF format, and FFN (nucleotide) and FAA (amino acid) FASTA sequence files. The outputs are multiple files of different types, including .ffn, .fna, and .fsa. The .ffn file is the main focus.

### BLAST (Basic Local Alignment Search Tool)
This tool identifies sequence similarity to a given reference set and identify sequence homology. This program can be run locally or on the NCBI website. The input is the 16S sequence that was determined from the PROKKA results. The output is the top hits in the database with the best taxonomic matches at the top.

### BWA mem
This tool alignis short reads to a reference sequence. The input to the program is a referece assembly and forward and reverse reads to map. The output is a SAM file that can be used to calculate coverage.

### Blobtools
THis tool visualizes the genome assembly and filters read and assembly data sets. The inputs are the contig fasta file, a "hits" file generated from the BLAST, and the SAM file generrated from the BWA program. The output includes a blobplot, which plots the GC, coverage, taxonomy, and contigs lengths on a single graph.

</details></details>

<details> <summary><H2> Results </H2></summary>

With the analysis above, the following can be performed:

<details> <summary><H3> Blast Hits </H3></summary>

![](https://github.com/ael1083/Toad-Protist/blob/main/images/BLAST%20Distribution.png?raw=true)

This shows the distribution of the top 30 BLAST hits on the 30 subject sequences against the query sequence. Most of the hits align at about 500bp, with a few hits also aligning at about 100bp. The top sequence (shown in purple), aligns at 100bp and then at 400-400bp, making that hit the closest match to the sample.

</details>

<details> <summary><H3> Genome Visualization </H3></summary>

![](https://github.com/user-attachments/assets/710b0667-e61c-42e7-99dc-f9329de3d574)

This graph consists of three separate parts. The top portion of the graph plots the GC proportion against the Span in kilobases. It can be seen that the Bufo have a GC proportion of about 0.5 and span up to about 1000kb. The bottom left portion plots the GC proportion against the coverage. It is seen that the coverage hovers at about 10<sup>-1</sup>. The bottom right portion of the graph plots the Span in kilobases against the coverage, and also hovers at around 10<sup>-1</sup>.

</details>

These results were not expected. It is assumed that a mistake was made during the creation of the SAM file.

</details></details>

<details> <summary><H2> Code </H2></summary>
  
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

#Index reference genome
bwa index contigs.fasta

#Map reads to construct SAM file
bwa mem -t 24 contigs.fasta 16S_sequences.fasta > raw_mapped.sam
less -S raw_mapped.sam

#Construct a coverage table
samtools view -@ 24 -Sb  raw_mapped.sam  | samtools sort -@ 24 -o sorted_mapp
ed.bam
samtools flagstat sorted_mapped.bam

#Index BAM file
samtools index sorted_mapped.bam

#Calculate per base coverage
bedtools genomecov -ibam sorted_mapped.bam > coverage.out
gen_input_table.py  --isbedfiles contigs.fasta coverage.out >  coverage_table
.tsv

#Create lookup table
blobtools create -i contigs.fasta -b sorted_mapped.bam -t contigs.fas
ta.vs.nt.cul5.maxt10.1e5.megablast.out -o blob_out

#Create output table & plot
blobtools view -i blob_out.blobDB.json -r all -o blob_taxonomy
grep -v '##' blob_taxonomy.blob_out.blobDB.table.txt
blobtools plot -i blob_out.blobDB.json -r genus

#Filter by length (<500bp) & coverage (>0)
mkdir mdibl-t3-2018-WGS/filtered_assembly
cp ../blob_taxonomy.blob_out.blobDB.table.txt ./
grep -v '#' blob_taxonomy.blob_out.blobDB.table.txt | awk -F'\t' '$2
< 5000' | awk -F'\t' '$5 > 0' | awk -F'\t' '{print $1}' > list_of_con
tigs_to_keep_len500_cov20.txt
less -S list_of_contigs_to_keep_len500_cov20.txt

#Filter assembly based on contigs list
filter_contigs_by_list.py ~/toad-analysis/trimmed-reads/spades_assemb
ly_default/contigs.fasta list_of_contigs_to_keep_len500_cov20.txt pro
tist_filtered.fasta

#Find average coverage
grep -f list_of_contigs_to_keep_len500_cov20.txt blob_taxonomy.blob_o
ut.blobDB.table.txt | awk '{w = w + $2; e = e + $5 * $2;} END {print
e/w}'

#BLAST final contigs
wget "https://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec"
blastn -reward 1 -penalty -5 -gapopen 3 -gapextend 3 -dust yes -soft_
masking true -evalue 700 -searchsp 1750000000000 -query protist_filte
red.fasta -subject UniVec  -outfmt 6 -out genome_vs_univec.6
```

</details></details>

<details> <summary><H2> Bibliography </H2></summary>

- This project was completed following the [MDIBL-T3-WGS-Tutorial](https://github.com/Joseph7e/MDIBL-T3-WGS-Tutorial?tab=readme-ov-file#organism-identification)
- [chatGPT](https://chatgpt.com/) helps explain what all the inputs do, outputs mean, and what can be done wtih them
- [NIH BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=MegaBlast&PROGRAM=blastn&PAGE_TYPE=BlastSearch&BLAST_SPEC=) was used to visualize the 16S sequences FASTA 

</details></details>
