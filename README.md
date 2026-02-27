# COMP383_PipelineProject 

Tools needed: Snakefile, Python, Biopython downloaded in python


First, I download the the samples from the two patients donors, at 2 days and 6 days post infection, using SRA numbers and converted to paired-end fastq files using wget, and store them to the folder full_data. 
code:
wget {link from SRA, normalization}
fasterqdump

I then stored test data in the folder test_data in this repo that will have the first 10,000 lines of the fastq files to make it easier to run them for the questions. 
code used for that: 
head -n 10000 SRR5660030_1.fastq > {path to file}


I then created the Snakefile and ran it with all nessary data downloaded.

I have stored longest_contig and question4.py in this repo, two python files that will be called in the snake file in the rules longest_contig and count_contig

To run snakefile on terminal:
snakemake --cores 4
