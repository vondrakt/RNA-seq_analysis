# RNA-seq_analysis
# Ref. based pipeline
# readscorrection of sequencing (optional)
rnaspades.py -s reads --only-error-correction -t 5 -o 

# trimming

wget http://genomedata.org/rnaseq-tutorial/illumina_multiplex.fa

trimmomatic PE -threads 30 -phred33 [corrected_input1].fastq.gz [corrected_input2].fastq.gz [trimmed_ourput]_1P.fastq.gz [trimmed_ourput]_1U.fastq.gz [trimmed_ourput]_2P.fastq.gz [trimmed_ourput]_2U.fastq.gz  ILLUMINACLIP:illumina_multiplex.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# fastqc and multiqc (optional)
fastqc [trimmed_ourput]_1P.fastq.gz [trimmed_ourput]_2P.fastq.gz
multiqc .

# STAR alignment 
wget https://www.arabidopsis.org/download_files/Genes/Araport11_genome_release/Araport11_GTF_genes_transposons.current.gtf.gz && gunzip Araport11_GTF_genes_transposons.current.gtf.gz
wget https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas.gz && gunzip  TAIR10_chr_all.fas.gz
STAR --runMode genomeGenerate --runThreadN 30 --genomeDir TAIR10_STAR --sjdbGTFfile Araport11_GTF_genes_transposons.current.gtf --genomeFastaFiles TAIR10_chr_all.fas
STAR --genomeDir TAIR10_STAR --runThreadN 30 --readFilesIn [trimmed_ourput]_1P.fastq.gz [trimmed_ourput]_2P.fastq.gz --readFilesCommand zcat --quantMode GeneCounts --outFileNamePrefix [output_name].sorted.bam --outSAMtype BAM SortedByCoordinate


# count reads and expression matrix based star BAM
# Create a stringtie samples file

sample 1	path_to_sample1.gtf

sample 2	path_to_sample2.gtf

sample 3	path_to_sample2.gtf



stringtie [output_name].sorted.bam -p 30 -G Araport11_GTF_genes_transposons.current.gtf -e -o [sample].gtf -A [sample].gene_abundances.tsv

wget https://ccb.jhu.edu/software/stringtie/dl/prepDE.py3

python3 prepDE.py3 -i sample_f.txt

# DE Analyis
# Create a samples.txt
condition1	condition1_R1
condition1	condition1_R2
condition1  condition1_R3  
condition2  condition2_R1
condition2  condition2_R2
condition2  condition2_R3




# edgeR, DESeq2 using Tinity DE_analysis script

run_DE_analysis.pl --matrix stringtie.gene/transcipt.counts.matrix --samples_file samples.txt --reference_sample [condition1/condition2] --method edgeR --output edgeR_genes
run_DE_analysis.pl --matrix stringtie.gene/transcipt.counts.matrix --samples_file samples.txt --reference_sample [condition1/condition2] --method DESeq2 --output DESeq2_genes

# OR
# R code from Protocol: Using StringTie with DESeq2
https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual 

# Extracting differentially expressed transcripts and generating heatmaps
# Extract those differentially expressed (DE) transcripts that are at least 4-fold differentially expressed at a significance of <= 0.05 in any of the pairwise sample comparisons:

analyze_diff_expr.pl --matrix stringtie.gene/transcipt.counts.matrix --samples samples.txt -P 0.05 -C 2 


# Functional Enrichment Analysis

Plaza (https://bioinformatics.psb.ugent.be/plaza/) 

TAIR GO Term Enrichment for Plants

PANTHER (http://pantherdb.org/webservices/go/overrep.jsp)
