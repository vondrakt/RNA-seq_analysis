#!/usr/bin/env python3
from optparse import OptionParser
import os
import re

# THIS IS A MASTER PYTHON SCRIPT FOR RUNNING A REFERENCE BASED RNA-SEQ ANALYSIS

# defining the arguments that need to be passed to the scripts
arguments = OptionParser()

arguments.add_option('-p', '--path', dest='path', help='path to wroking directory')
arguments.add_option('-s', '--samples', dest='samples', help='sample IDs, separated by comma')
arguments.add_option('-r', '--replicates', dest='replicates', help='replicate IDs, separated by comma')
arguments.add_option('-g', '--genome', dest='genome', help='reference genome to be used for the analysis')
arguments.add_option('-a', '--annotation-gtf', dest='annotation_gtf', help='annotation of the reference genome in gtf format')
arguments.add_option('-b', '--annotation-gff', dest='annotation_gff', help='annotation of the reference genome in gff format')
arguments.add_option('-c', '--config-file', dest='config_DE', help='configuration file for the DE analysis')

(options, args) = arguments.parse_args()
if options.path is None or options.samples is None or options.replicates is None\
        or options.genome is None or options.annotation_gtf is None or options.annotation_gff is None or options.config_DE is None:
    # if one of the arguments is missing
    print('\n----------> A mandatory option is missing !\n')  # raise an error
    arguments.print_help()  # and provide the help
    exit(-1)  # exit the script
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Due to incompatibility with rnaspades, this step is skipped

# Processing names
samples = options.samples.split(',')
replicates = options.replicates.split(',')

for replicate in replicates:
    name1 = options.path+'/'+samples[0]+replicate
    name2 = options.path+'/'+samples[1]+replicate
    print(name1, name2)
    # Trimming the reads
    trimming_command = 'trimmomatic PE -threads 30 -phred33 %s.fastq.gz' \
                       ' %s.fastq.gz %s_1P.fastq.gz %s_1U.fastq.gz %s_2P.fastq.gz %s_2U.fastq.gz ' \
                       'ILLUMINACLIP:illumina_multiplex.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36' % (name1, name2, name1, name1, name2, name2)
    print(trimming_command)
    os.system(trimming_command)

    # Fastqc
    fastqc_command = 'fastqc %s_1P.fastq.gz %s_2P.fastq.gz' % (name1, name2)
    multiqc_command = 'multiqc %s_1P.fastq.gz %s_2P.fastq.gz' % (name1, name2)
    print(fastqc_command)
    print(multiqc_command)
    os.system(fastqc_command)
    os.system(multiqc_command)

    # STAR alignment
    STAR_command_reference = 'STAR --runMode genomeGenerate --runThreadN 30 --genomeDir TAIR10_STAR --sjdbGTFfile %s --genomeFastaFiles %s' % (options.annotation_gtf, options.genome)
    print(STAR_command_reference)
    os.system(STAR_command_reference)
    STAR_command_reads_sample1 = 'STAR --genomeDir TAIR10_STAR --runThreadN 10 --readFilesIn %s_1P.fastq.gz ' \
                         '--readFilesCommand zcat --quantMode GeneCounts --outFileNamePrefix %s.sorted.bam -' \
                         '-outSAMtype BAM SortedByCoordinate' % (samples[0]+replicate, samples[0]+replicate)

    STAR_command_reads_sample2 = 'STAR --genomeDir TAIR10_STAR --runThreadN 10 --readFilesIn %s_2P.fastq.gz ' \
                         '--readFilesCommand zcat --quantMode GeneCounts --outFileNamePrefix %s.sorted.bam -' \
                         '-outSAMtype BAM SortedByCoordinate' % (samples[1] + replicate, samples[1] + replicate)

    print(STAR_command_reads_sample1)
    print(STAR_command_reads_sample2)
    os.system(STAR_command_reads_sample1)
    os.system(STAR_command_reads_sample2)

    # count reads and expression matrix based star BAM
    # Create a stringtie samples file
    stringtie_command_sample1 = 'stringtie %s.sorted.bamAligned.sortedByCoord.out.bam -p 30 -G %s -e -o %s.gtf -A %s.gene_abundances.tsv' % \
                        (samples[0]+replicate, options.annotation_gff, samples[0]+replicate, samples[0]+replicate,)

    stringtie_command_sample2 = 'stringtie %s.sorted.bamAligned.sortedByCoord.out.bam -p 30 -G %s -e -o %s.gtf -A %s.gene_abundances.tsv' % \
                                (samples[1] + replicate, options.annotation_gff, samples[1] + replicate,
                                 samples[1] + replicate,)
    print(stringtie_command_sample1)
    print(stringtie_command_sample2)
    os.system(stringtie_command_sample1)
    os.system(stringtie_command_sample2)

# Create the samples.txt file
samples_txt = open('prepDE_samples.txt', 'w')

for sample in samples:
    for replicate in replicates:
        list = [sample+replicate, options.path+'/'+sample+replicate+'.gtf']
        list = '\t'.join(list)
        samples_txt.write(list+'\n')
        #print(list)
        #count += 1

samples_txt.close()
# Generating count matrices
command = 'prepDE.py -i prepDE_samples.txt'
print(command)
os.system(command)

# DE analysis
# The configuration file for the DE analysis must be in format
# condition1    condition1_R1
# condition1    condition1_R2
# condition1    condition1_R3
# condition2    condition2_R1
# condition2    condition2_R2
# condition2    condition2_R3

# first convert the .csv file to a .tsv file
#convert_command = 'cat gene_count_matrix.csv | perl -pe "s/,/\t/g > gene_count_matrix.tsv"'
out_tsv = open('./gene_count_matrix.tsv', 'w')
with open('./gene_count_matrix.csv') as c:
    for line in c:
        line = re.sub(',', '\t', line)
        out_tsv.write(line)
out_tsv.close()

DE_analysis_command_edgeR = 'run_DE_analysis.pl --matrix gene_count_matrix.tsv --samples_file %s ' \
                      '--reference_sample condition1 --method edgeR --output edgeR_genes' % options.config_DE
os.system(DE_analysis_command_edgeR)

DE_analysis_command_DESeq2 = 'run_DE_analysis.pl --matrix gene_count_matrix.tsv --samples_file %s ' \
                      '--reference_sample condition1 --method DESeq2 --output DESeq2_genes' % options.config_DE
os.system(DE_analysis_command_DESeq2)

# Extracting differentially expressed transcripts and generating heatmaps
# Extract those differentially expressed (DE) transcripts that are at least 4-fold differentially expressed
# at a significance of <= 0.05 in any of the pairwise sample comparisons:

analyze_DE = 'analyze_diff_expr.pl --matrix transcript_count_matrix.csv --samples samples.txt -P 0.05 -C 2'
os.system(analyze_DE)











