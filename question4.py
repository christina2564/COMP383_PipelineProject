import sys
import argparse
from Bio import SeqIO # uses seqIO to parse thorugh 
#function to parse command line arguments
def check_arg(args=None):
    parser = argparse.ArgumentParser(description="<add description of what script does>")
    parser.add_argument("-i", "--input", 
    help="input file", 
    required=True)
    parser.add_argument("-o", "--output",
    help="output file",
    required=True)
    parser.add_argument("-s", "--sample", # adds third reuired argunt to print sample name in report 
    help="sample name",
    required=True)
    return parser.parse_args(args)

#retrieve command line arguments
arguments = check_arg(sys.argv[1:])
infile = arguments.input
outfile = arguments.output
sample_name = arguments.sample # gives you sample id

num_contigs = 0 # intialize number of contigs to 0
total_bp = 0 # intitalizes total bp to 0

for f in SeqIO.parse(infile,"fasta"): # for each sequence record in the fasta file it will parse thoruh, and loops through file 
    if len(f.seq) > 1000: # if the length of the sequence in bigger than 1000
        num_contigs+= 1 # then it adds to the contig_count by 1 
        total_bp += len(f.seq) # adds length of sequences to total_bp

with open(outfile, "w") as out: # writing the output file 
    out.write(f"In the assembly of sample {sample_name}, there are {num_contigs} contgis > 1000 bp and {total_bp} total bp. \n") # formats the report 