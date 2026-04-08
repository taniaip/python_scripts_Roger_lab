------------FOR TRIMMAMATIC:-----------
#!/bin/bash

#$ -S /bin/bash
#$ -cwd
#$ -pe threaded 8
#$ -q 256G-batch,768G-batch,2T-batch

THREADS=8

# Input raw paired-end RNA-seq reads
R1=/path/to/raw_forward_reads.fastq.gz
R2=/path/to/raw_reverse_reads.fastq.gz

# Output trimmed paired and unpaired reads
FW_PRD=trimmed_forward_paired.fastq.gz
FW_UNP=trimmed_forward_unpaired.fastq.gz
RV_PRD=trimmed_reverse_paired.fastq.gz
RV_UNP=trimmed_reverse_unpaired.fastq.gz

source activate trimmomatic

trimmomatic PE -threads $THREADS \
    $R1 $R2 \
    $FW_PRD $FW_UNP $RV_PRD $RV_UNP \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:1:true \
    LEADING:3 \
    TRAILING:3 \
    MINLEN:60 \
    AVGQUAL:30

conda deactivate


------------FOR HISAT2:------------
#!/bin/bash

#$ -S /bin/bash
. /etc/profile
#$ -cwd
#$ -pe threaded 10

# Directory containing input/output files
base='/path/to/working_directory'

# Input trimmed RNA-seq reads
FW_PRD=trimmed_forward_paired.fastq.gz
FW_UNP=trimmed_forward_unpaired.fastq.gz
RV_PRD=trimmed_reverse_paired.fastq.gz
RV_UNP=trimmed_reverse_unpaired.fastq.gz

# Input reference genome assembly
ASSEMBLY='reference_genome.fasta'

# Output alignment prefix and HISAT2 index basename
PREFIX=rnaseq_alignment_output
INDEX=reference_genome_index

# Load HISAT2
source activate hisat2

# Build genome index
hisat2-build -f $ASSEMBLY $INDEX

# Align RNA-seq reads to the genome
# PAY ATTENTION to choose one of RF/FR!!
hisat2 \
    -q --threads 10 --phred33 \
    --max-intronlen 1000 \
    --rna-strandness FR/RF -k 2 \
    -x $INDEX -1 $FW_PRD -2 $RV_PRD -U $FW_UNP,$RV_UNP | \
    samtools sort --threads 10 -o $PREFIX.sort.bam -

# Index BAM file
samtools index -@ 10 $PREFIX.sort.bam $PREFIX.sort.bam.bai

conda deactivate
