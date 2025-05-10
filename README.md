# Toad-Protist Project

Alexandria Lyden, Lawrence Gordon, W. Kelley Thomas

<details> <summary><H2> Background </H2></summary>

The data for this analysis was provided by the Hubbard Center for Genome Studies at the University of New Hampshire. It consists of 16s data in paired-end 250 bp reads that were amplified by Illumina HiSeq 2500. The files were made up of 1 sample taken from the Anaxyrus americanus, known as the American Toad. With this data, the goal was to assemble the genome of the protist, identify the protist's class, and assess the genome. The protist, Amphibiothecum penneri, was identified from [this](https://pubmed.ncbi.nlm.nih.gov/16456158/) paper.

</details></details>

<details> <summary><H2> Methods </H2></summary>

### Trimmomatic

### SPADes

### PROKKA

### QUAST

### BUSCO

### BLAST

### BWA mem

### Blobtools

</details></details>

<details> <summary><H2> Results </H2></summary>

With the analysis above, the following can be performed:

<details> <summary><H3> Blast Hits </H3></summary>

![](https://github.com/ael1083/Toad-Protist/blob/main/images/BLAST%20Distribution.png?raw=true)

Explain

</details>

<details> <summary><H3> Genome Visualization </H3></summary>

![](https://github.com/user-attachments/assets/710b0667-e61c-42e7-99dc-f9329de3d574)

Explain

</details>

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
