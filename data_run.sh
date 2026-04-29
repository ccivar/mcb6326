#!/bin/bash
#SBATCH --job-name=data_run # Job name
#SBATCH --mail-type=END,FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user= # Where to send mail
#SBATCH --ntasks=16 # Run on 16 single CPUs
#SBATCH --mem=6gb # Job memory request
#SBATCH --time=24:00:00 # Time limit hrs:min:sec
#SBATCH --output=serial_test_%j.log # Standard output and error log

#load modules needed
ml sra
ml fastqc
ml trimmomatic 
ml hisat2
ml samtools
ml htseq

####### Get Resistant Control Samples #########

#get data Brassica Resistant control Sample 1
wget  -O Brassica_Resistant_control_1 https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR14730443/SRR14730443

#split SRA file into 2 reads then remove file
fastq-dump --split-files Brassica_Resistant_control_1
rm Brassica_Resistant_control_1

#zip the fastq files
gzip Brassica_Resistant_control_1_1.fastq
gzip Brassica_Resistant_control_1_2.fastq


#get data Brassica Resistant control Sample 2
wget  -O Brassica_Resistant_control_2 https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR14730442/SRR14730442

#split SRA file into 2 reads then remove file
fastq-dump --split-files Brassica_Resistant_control_2
rm Brassica_Resistant_control_2

#zip the fastq files
gzip Brassica_Resistant_control_2_1.fastq
gzip Brassica_Resistant_control_2_2.fastq



####### Get Resistant Infected Samples #########

#get data Brassica Resistant Infected Sample 1
wget  -O Brassica_Resistant_infected_1 https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR14730433/SRR14730433

#split SRA file into 2 reads then remove file
fastq-dump --split-files Brassica_Resistant_infected_1
rm Brassica_Resistant_infected_1

#zip the fastq files
gzip Brassica_Resistant_infected_1_1.fastq
gzip Brassica_Resistant_infected_1_2.fastq


#get data Brassica Resistant Infected Sample 2
wget  -O Brassica_Resistant_infected_2 https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR14730432/SRR14730432

#split SRA file into 2 reads then remove file
fastq-dump --split-files Brassica_Resistant_infected_2
rm Brassica_Resistant_infected_2

#zip the fastq files
gzip Brassica_Resistant_infected_2_1.fastq
gzip Brassica_Resistant_infected_2_2.fastq





####### Get Susceptible Control Samples #########

#get data Brassica Susceptible control Sample 1
wget  -O Brassica_Susceptible_control_1 https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR14730426/SRR14730426

#split SRA file into 2 reads then remove file
fastq-dump --split-files Brassica_Susceptible_control_1
rm Brassica_Susceptible_control_1

#zip the fastq files
gzip Brassica_Susceptible_control_1_1.fastq
gzip Brassica_Susceptible_control_1_2.fastq


#get data Brassica Susceptible control Sample 2
wget  -O Brassica_Susceptible_control_2 https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR14730425/SRR14730425

#split SRA file into 2 reads then remove file
fastq-dump --split-files Brassica_Susceptible_control_2
rm Brassica_Susceptible_control_2

#zip the fastq files
gzip Brassica_Susceptible_control_2_1.fastq
gzip Brassica_Susceptible_control_2_2.fastq


####### Get Infected Susceptible Samples #########

#get data Brassica Susceptible Infected Sample 1
wget  -O Brassica_Susceptible_infected_1 https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR14730429/SRR14730429

#split SRA file into 2 reads then remove file
fastq-dump --split-files Brassica_Susceptible_infected_1
rm Brassica_Susceptible_infected_1

#zip the fastq files
gzip Brassica_Susceptible_infected_1_1.fastq
gzip Brassica_Susceptible_infected_1_2.fastq


#get data Brassica Susceptible Infected Sample 2
wget  -O Brassica_Susceptible_infected_2 https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR14730428/SRR14730428

#split SRA file into 2 reads then remove file
fastq-dump --split-files Brassica_Susceptible_infected_2
rm Brassica_Susceptible_infected_2

#zip the fastq files
gzip Brassica_Susceptible_infected_2_1.fastq
gzip Brassica_Susceptible_infected_2_2.fastq


#check quality of the reads
#fastqc -t 16 *.fastq.gz


#build index 
hisat2-build -f Brassica_Genome_DNA.fa Brassica_Genome


#Resistant Controls
trimmomatic PE -phred33 -threads 16 Brassica_Resistant_control_1_1.fastq.gz \
Brassica_Resistant_control_1_2.fastq.gz \
Brassica_Resistant_control_1_1_paired_trimmed.fq.gz Brassica_Resistant_control_1_1_unpaired_trimmed.fq.gz \
Brassica_Resistant_control_1_2_paired_trimmed.fq.gz Brassica_Resistant_control_1_2_unpaired_trimmed.fq.gz \
ILLUMINACLIP:${HPC_TRIMMOMATIC_ADAPTER}/TruSeq3-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:75

trimmomatic PE -phred33 -threads 16 Brassica_Resistant_control_2_1.fastq.gz \
Brassica_Resistant_control_2_2.fastq.gz \
Brassica_Resistant_control_2_1_paired_trimmed.fq.gz Brassica_Resistant_control_2_1_unpaired_trimmed.fq.gz \
Brassica_Resistant_control_2_2_paired_trimmed.fq.gz Brassica_Resistant_control_2_2_unpaired_trimmed.fq.gz \
ILLUMINACLIP:${HPC_TRIMMOMATIC_ADAPTER}/TruSeq3-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:75


#Susceptible Controls
trimmomatic PE -phred33 -threads 16 Brassica_Susceptible_control_1_1.fastq.gz \
Brassica_Susceptible_control_1_2.fastq.gz \
Brassica_Susceptible_control_1_1_paired_trimmed.fq.gz Brassica_Susceptible_control_1_1_unpaired_trimmed.fq.gz \
Brassica_Susceptible_control_1_2_paired_trimmed.fq.gz Brassica_Susceptible_control_1_2_unpaired_trimmed.fq.gz \
ILLUMINACLIP:${HPC_TRIMMOMATIC_ADAPTER}/TruSeq3-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:75

trimmomatic PE -phred33 -threads 16 Brassica_Susceptible_control_2_1.fastq.gz \
Brassica_Susceptible_control_2_2.fastq.gz \
Brassica_Susceptible_control_2_1_paired_trimmed.fq.gz Brassica_Susceptible_control_2_1_unpaired_trimmed.fq.gz \
Brassica_Susceptible_control_2_2_paired_trimmed.fq.gz Brassica_Susceptible_control_2_2_unpaired_trimmed.fq.gz \
ILLUMINACLIP:${HPC_TRIMMOMATIC_ADAPTER}/TruSeq3-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:75


#Resistant Infected
trimmomatic PE -phred33 -threads 16 Brassica_Resistant_infected_1_1.fastq.gz \
Brassica_Resistant_infected_1_2.fastq.gz \
Brassica_Resistant_infected_1_1_paired_trimmed.fq.gz Brassica_Resistant_infected_1_1_unpaired_trimmed.fq.gz \
Brassica_Resistant_infected_1_2_paired_trimmed.fq.gz Brassica_Resistant_infected_1_2_unpaired_trimmed.fq.gz \
ILLUMINACLIP:${HPC_TRIMMOMATIC_ADAPTER}/TruSeq3-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:75

trimmomatic PE -phred33 -threads 16 Brassica_Resistant_infected_2_1.fastq.gz \
Brassica_Resistant_infected_2_2.fastq.gz \
Brassica_Resistant_infected_2_1_paired_trimmed.fq.gz Brassica_Resistant_infected_2_1_unpaired_trimmed.fq.gz \
Brassica_Resistant_infected_2_2_paired_trimmed.fq.gz Brassica_Resistant_infected_2_2_unpaired_trimmed.fq.gz \
ILLUMINACLIP:${HPC_TRIMMOMATIC_ADAPTER}/TruSeq3-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:75

#Susceptible Infected
trimmomatic PE -phred33 -threads 16 Brassica_Susceptible_infected_1_1.fastq.gz \
Brassica_Susceptible_infected_1_2.fastq.gz \
Brassica_Susceptible_infected_1_1_paired_trimmed.fq.gz Brassica_Susceptible_infected_1_1_unpaired_trimmed.fq.gz \
Brassica_Susceptible_infected_1_2_paired_trimmed.fq.gz Brassica_Susceptible_infected_1_2_unpaired_trimmed.fq.gz \
ILLUMINACLIP:${HPC_TRIMMOMATIC_ADAPTER}/TruSeq3-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:75

trimmomatic PE -phred33 -threads 16 Brassica_Susceptible_infected_2_1.fastq.gz \
Brassica_Susceptible_infected_2_2.fastq.gz \
Brassica_Susceptible_infected_2_1_paired_trimmed.fq.gz Brassica_Susceptible_infected_2_1_unpaired_trimmed.fq.gz \
Brassica_Susceptible_infected_2_2_paired_trimmed.fq.gz Brassica_Susceptible_infected_2_2_unpaired_trimmed.fq.gz \
ILLUMINACLIP:${HPC_TRIMMOMATIC_ADAPTER}/TruSeq3-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:75


#map reads to the Brassica Genome and convert to bam files

hisat2 -p 16 -x Brassica_Genome -1 Brassica_Resistant_control_1_1_paired_trimmed.fq.gz -2 Brassica_Resistant_control_1_2_paired_trimmed.fq.gz -S resistant_control_rep1.sam
samtools view -bS resistant_control_rep1.sam > resistant_control_rep1.bam
rm *.sam

hisat2 -p 16 -x Brassica_Genome -1 Brassica_Resistant_control_2_1_paired_trimmed.fq.gz -2 Brassica_Resistant_control_2_2_paired_trimmed.fq.gz -S resistant_control_rep2.sam
samtools view -bS resistant_control_rep2.sam > resistant_control_rep2.bam
rm *.sam

hisat2 -p 16 -x Brassica_Genome -1 Brassica_Susceptible_control_1_1_paired_trimmed.fq.gz -2 Brassica_Susceptible_control_1_2_paired_trimmed.fq.gz -S susceptible_control_rep1.sam
samtools view -bS susceptible_control_rep1.sam > susceptible_control_rep1.bam
rm *.sam

hisat2 -p 16 -x Brassica_Genome -1 Brassica_Susceptible_control_2_1_paired_trimmed.fq.gz -2 Brassica_Susceptible_control_2_2_paired_trimmed.fq.gz -S susceptible_control_rep2.sam
samtools view -bS susceptible_control_rep2.sam > susceptible_control_rep2.bam
rm *.sam

hisat2 -p 16 -x Brassica_Genome -1 Brassica_Resistant_infected_1_1_paired_trimmed.fq.gz -2 Brassica_Resistant_infected_1_2_paired_trimmed.fq.gz -S resistant_infected_rep1.sam
samtools view -bS resistant_infected_rep1.sam > resistant_infected_rep1.bam
rm *.sam

hisat2 -p 16 -x Brassica_Genome -1 Brassica_Resistant_infected_2_1_paired_trimmed.fq.gz -2 Brassica_Resistant_infected_2_2_paired_trimmed.fq.gz -S resistant_infected_rep2.sam
samtools view -bS resistant_infected_rep2.sam > resistant_infected_rep2.bam
rm *.sam

hisat2 -p 16 -x Brassica_Genome -1 Brassica_Susceptible_infected_1_1_paired_trimmed.fq.gz -2 Brassica_Susceptible_infected_1_2_paired_trimmed.fq.gz -S susceptible_infected_rep1.sam
samtools view -bS susceptible_infected_rep1.sam > susceptible_infected_rep1.bam
rm *.sam

hisat2 -p 16 -x Brassica_Genome -1 Brassica_Susceptible_infected_2_1_paired_trimmed.fq.gz -2 Brassica_Susceptible_infected_2_2_paired_trimmed.fq.gz -S susceptible_infected_rep2.sam
samtools view -bS susceptible_infected_rep2.sam > susceptible_infected_rep2.bam
rm *.sam


#sort by coordinate
samtools sort -o resistant_control_rep1_sorted.bam resistant_control_rep1.bam
samtools sort -o resistant_control_rep2_sorted.bam resistant_control_rep2.bam
samtools sort -o susceptible_control_rep1_sorted.bam susceptible_control_rep1.bam
samtools sort -o susceptible_control_rep2_sorted.bam susceptible_control_rep2.bam
samtools sort -o resistant_infected_rep1_sorted.bam resistant_infected_rep1.bam
samtools sort -o resistant_infected_rep2_sorted.bam resistant_infected_rep2.bam
samtools sort -o susceptible_infected_rep1_sorted.bam susceptible_infected_rep1.bam
samtools sort -o susceptible_infected_rep2_sorted.bam susceptible_infected_rep2.bam

rm  resistant_control_rep1.bam
rm resistant_control_rep2.bam
rm susceptible_control_rep1.bam
rm susceptible_control_rep2.bam
rm resistant_infected_rep1.bam
rm resistant_infected_rep2.bam
rm susceptible_infected_rep1.bam
rm susceptible_infected_rep2.bam

#index
samtools index resistant_control_rep1_sorted.bam
samtools index resistant_control_rep2_sorted.bam
samtools index susceptible_control_rep1_sorted.bam
samtools index susceptible_control_rep2_sorted.bam
samtools index resistant_infected_rep1_sorted.bam
samtools index resistant_infected_rep2_sorted.bam
samtools index susceptible_infected_rep1_sorted.bam
samtools index susceptible_infected_rep2_sorted.bam


#count reads

#htseq-count -f bam -r pos -m intersection-nonempty -s no -t gene -i ID -o resistant_control_rep1_counts resistant_control_rep1_sorted.bam  Brassica_oleracea.genes.gff3 > resistant_control_rep1_gene_summary.txt
rm resistant_control_rep1_counts
htseq-count -f bam -r pos -m intersection-nonempty -s no -t gene -i ID -o resistant_control_rep2_counts resistant_control_rep2_sorted.bam  Brassica_oleracea.genes.gff3 > resistant_control_rep2_gene_summary.txt
rm resistant_control_rep2_counts
htseq-count -f bam -r pos -m intersection-nonempty -s no -t gene -i ID -o susceptible_control_rep1_counts susceptible_control_rep1_sorted.bam  Brassica_oleracea.genes.gff3 > susceptible_control_rep1_gene_summary.txt
rm susceptible_control_rep1_counts
htseq-count -f bam -r pos -m intersection-nonempty -s no -t gene -i ID -o susceptible_control_rep2_counts susceptible_control_rep2_sorted.bam  Brassica_oleracea.genes.gff3 > susceptible_control_rep2_gene_summary.txt
rm susceptible_control_rep2_counts

htseq-count -f bam -r pos -m intersection-nonempty -s no -t gene -i ID -o resistant_infected_rep1_counts resistant_infected_rep1_sorted.bam  Brassica_oleracea.genes.gff3 > resistant_infected_rep1_gene_summary.txt
rm resistant_infected_rep1_counts
htseq-count -f bam -r pos -m intersection-nonempty -s no -t gene -i ID -o resistant_infected_rep2_counts resistant_infected_rep2_sorted.bam  Brassica_oleracea.genes.gff3 > resistant_infected_rep2_gene_summary.txt
rm resistant_infected_rep2_counts
htseq-count -f bam -r pos -m intersection-nonempty -s no -t gene -i ID -o susceptible_infected_rep1_counts susceptible_infected_rep1_sorted.bam  Brassica_oleracea.genes.gff3 > susceptible_infected_rep1_gene_summary.txt
rm susceptible_infected_rep1_counts
htseq-count -f bam -r pos -m intersection-nonempty -s no -t gene -i ID -o susceptible_infected_rep2_counts susceptible_infected_rep2_sorted.bam  Brassica_oleracea.genes.gff3 > susceptible_infected_rep2_gene_summary.txt
rm susceptible_infected_rep2_counts
#remove files

rm *.sam
rm *.bam
rm *fq.gz
rm *.bai