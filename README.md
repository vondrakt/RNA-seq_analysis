# Running the reference based RNA-seq analysis

Create a conda enviroment and install necessary dependecies

`mamba env create -n RNA-seq_analysis -f DE_analysis.yaml`

Activate the conda enviroment

`conda activate RNA-seq_analysis`

Running the master python script

`./master_RNA_seq.py  -p /path/to/directory -s 17NQ004_P1_R11,17NQ004_P2_R5 -r A,B,C -g TAIR10_chr_all.fas -a Araport11_GTF_genes_transposons.current.gtf -b TAIR10_GFF3_genes_transposons.gff -c config_DE.txt`

The files used for the commands are in test_data.zip.
The reference Arabidopsis genome, as well as the .gff and .gtf annotation must be downloaded from [https://www.arabidopsis.org/]

Choosing between analyzing the gene or transcript counts can be done with the master_RNA_seq_02.py script.
The option -x is used to determine which of the two options will be used.

`./master_RNA_seq_02.py -p /path/to/directory -s 17NQ004_P1_R11,17NQ004_P2_R5 -r A,B,C -g TAIR10_chr_all.fas -a Araport11_GTF_genes_transposons.current.gtf -b TAIR10_GFF3_genes_transposons.gff -c config_DE.txt -x TRANSCRIPT|GENE`
