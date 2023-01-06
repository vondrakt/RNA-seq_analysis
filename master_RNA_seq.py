#!/usr/bin/env python3
from optparse import OptionParser
import os

# THIS IS A MASTER PYTHON SCRIPT FOR RUNNING A REFERENCE BASED RNA-SEQ ANALYSIS

# defining the arguments that need to be passed to the scripts
arguments = OptionParser()

arguments.add_option('-p', '--path', dest='path', help='path to wroking directory')
arguments.add_option('-s', '--samples', dest='samples', help='sample IDs, separated by comma')
arguments.add_option('-r', '--replicates', dest='replicates', help='replicate IDs, separated by comma')
arguments.add_option('-g', '--genome', dest='genome', help='reference genome to be used for the analysis')
arguments.add_option('-a', '--annotation', dest='annotation', help='annotation of the reference genome')
arguments.add_option('-o', '--out', dest='out', help='path to output file')

(options, args) = arguments.parse_args()
if options.path is None or options.samples is None or options.replicates is None\
        or options.out is None or options.genome is None or options.annotation is None:
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
    #os.system(trimming_command)

    # Fastqc
    fastqc_command = 'fastqc %s_1P.fastq.gz %s_2P.fastq.gz' % (name1, name2)
    multiqc_command = 'multiqc %s_1P.fastq.gz %s_2P.fastq.gz' % (name1, name2)
    print(fastqc_command)
    print(multiqc_command)
    # os.system(fastqc_command)
    # os.system(multiqc_command)

    # STAR alignment
    STAR_command_reference = 'STAR --runMode genomeGenerate --runThreadN 30 --genomeDir TAIR10_STAR --sjdbGTFfile %s --genomeFastaFiles %s' % (options.annotation, options.genome)
    print(STAR_command_reference)
    # os.system(STAR_command)
    STAR_command_reads = 'STAR --genomeDir TAIR10_STAR --runThreadN 30 --readFilesIn %s_1P.fastq.gz %s_2P.fastq.gz ' \
                         '--readFilesCommand zcat --quantMode GeneCounts --outFileNamePrefix %s_vs_%s.sorted.bam -' \
                         '-outSAMtype BAM SortedByCoordinate' % (name1, name2, name1, name2)
    print(STAR_command_reads)
    # os.system(STAR_command_reads)

    # count reads and expression matrix based star BAM
    # Create a stringtie samples file




