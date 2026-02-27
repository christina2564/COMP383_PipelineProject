
samples = ["SRR5660030", "SRR5660033", "SRR5660044", "SRR5660045"] # make list called samples, storing 4 sample names, will be used in {sample}
rule all: # outputs used, all the ends, final targets used, the leaf nodes and some intermediates
    input:
        expand("results/{sample}.before.txt", sample= samples), # expand() makes many file names by replacing {sample} with the items in the list samples
        expand("results/{sample}.sam", sample=samples),
        expand("results/{sample}_mapped.1.fastq",sample=samples),
        expand("results/{sample}_mapped.2.fastq", sample=samples),
        expand("assemblies/{sample}/contigs.fasta", sample=samples),
        "Sebastian_PipelineReport.txt", # our final report 
        expand("results/{sample}.report_q2.txt", sample=samples),
        expand("results/{sample}_stats.txt", sample=samples),
        expand("results/{sample}_blast.txt", sample=samples)



rule retrieve_reference_genome: # rule to retrieve 
    output:
        "reference/hcmv.fna" # the final file we want 
    shell:
        """
        mkdir -p reference # makes reference dictionary, -p makes parent dictionaries if needed
        datasets download genome accession GCF_000845245.1 --include gff3,rna,cds,protein,genome,seq-report --filename reference_dataset.zip  # download the reference genome using accession number, renamed it because we use ncbi.dataset in other rule as well 
        unzip reference_dataset.zip -d reference_tmp  # unzips file and puts it dictonary called reference_tmp
        mv reference_tmp/ncbi_dataset/data/GCF_000845245.1/GCF_000845245.1_ViralProj14559_genomic.fna {output} # moves the fna file we wnat to output/ reference and renames if hcmv.fna to make it easier

        rm -r reference_tmp reference_dataset.zip # remove reference_tmp dictionary and the zip file
         
        """
rule create_index1: # we build an index from reference genome to use in bowtie2
    input: 
        "reference/hcmv.fna" # input is the file made from before rule, depends on before rule
    output: # listing all the index files we get from this rule
        "reference/hcmv.1.bt2", 
        "reference/hcmv.2.bt2", 
        "reference/hcmv.3.bt2", 
        "reference/hcmv.4.bt2",
        "reference/hcmv.rev.1.bt2",
        "reference/hcmv.rev.2.bt2"
        
    shell: 
        "bowtie2-build {input} reference/hcmv" # creates index files with prefix from reference/hcmv from input


rule count_reads_before: # counts how many readfs are in each sample before bowtie 
    input:
        R1="test_data/{sample}_1_test.fastq" # # name our input R1, and uses {sample} to call each of the samples
    output:
        "results/{sample}.before.txt" # txt file we get with number of reads, rule will make one for ecah sample
    shell:
        "echo $(( $(cat {input.R1} | wc -l) / 4 )) > {output}"  # we print out the entire file using cat, and then count each line using wc -l and then divide by 4 becuase fastq files have 4 lines per read, using $() lets us do the math, echo prints the result, > will write to the output file


 

rule bowtie2: # rule to align sample reads to reference genome 
    input: 
        R1="test_data/{sample}_1_test.fastq", # each sample's forward read file 
        R2="test_data/{sample}_2_test.fastq", # each sample's reverse read file 
        idx=rules.create_index1.output # tells snakmake this rule depends on outputs from create_index1, so brings us all the ouput from there, and will run rules in right order 
    output:
        sam ="results/{sample}.sam", # reults in sam format
        mapped1="results/{sample}_mapped.1.fastq", #forward reads that algined 
        mapped2="results/{sample}_mapped.2.fastq" # fastq with reverse reads tha taligned
    shell:
        "bowtie2 --quiet -x reference/hcmv -1 {input.R1} -2 {input.R2} --al-conc results/{wildcards.sample}_mapped.fastq -S {output.sam}" # -x to look for reverence/hcmv, -1 reads forward file, -2 reads reverse file, --al-conc writes algined paired reads to seperate files, so we get reads that mapped to index only, -S writes sam file to that file name 



rule count_reads_after:
    input:
        "results/{sample}_mapped.1.fastq" # results from bowtie -al-conc
    output:
        "results/{sample}.after.txt" # output with number of mapped reads
    shell: 
        "echo $(($(wc -l < {input})/ 4 )) > {output}"  # counts lines in fastq file and divde by 4 

rule q1_report: 
    input:
        before="results/{sample}.before.txt", # file from count_reads_before 
        after="results/{sample}.after.txt" # file from count_reads_after
    output:
        report="results/{sample}.report_q2.txt" # anwsers question 2   
    shell: 
        """
        before=$(cat {input.before}) # prints out content from the txt, and stores it in this variable 
        after=$(cat {input.after}) # prints out contntn froma fter file, and stores it 
        echo "Sample {wildcards.sample} had $before read pairs before  and $after read pairs after " > {output} # # writes out the summary of results and puts into output file 
        """

#3

rule spades_run: #runs SPAdes
    input: # mapped reads from bowtie 
        R1="results/{sample}_mapped.1.fastq",
        R2="results/{sample}_mapped.2.fastq"
    output: # contigs we get 
        "assemblies/{sample}/contigs.fasta"
    shell:
        "spades.py -k 99 -t 4 --only-assembler -1 {input.R1} -2 {input.R2} -o assemblies/{wildcards.sample}" # running spades at k-mer size 99, usings 4 threads, skips to assembly because we have mapped reads, reads forard and reverse files, -o sets outptu directory out from output, and we just want contigs.fasta 

#4
rule count_contigs: # counting contigs 
    input:
        "assemblies/{sample}/contigs.fasta" # using ouut from spades_run
    output:
        "results/{sample}_stats.txt" # get a files with the summary of the output 
    shell:
        "python question4.py -i {input} -o {output} -s {wildcards.sample}" # calls python file 

#5
rule longest_contig:
    input:
        "assemblies/{sample}/contigs.fasta" # also uses output from spades_run
    output: 
        "assemblies/{sample}/longest_contig.fna" # prints out the fasta file with the longest contig 
    shell:
        "python3 longest_contig -i {input} -o {output} " # runs the python file that finds us the longest contigs 

rule download_datasets:
    output:
        directory("bhv_dataset") # make new folder/dcitionary that will holds all genomes from ncbi Betaherpesvirinae
    shell:
        """
        datasets download virus genome taxon Betaherpesvirinae --refseq --include genome --filename bhv_dataset.zip # downloods the genome of the virus family, and only refseq genomes, only include genomes fasta file,and renames it from ncbi_dataset to make ti less confusing 
        unzip bhv_dataset.zip -d bhv_dataset # unzips files and extracts into folder bhv_dataset
        rm bhv_dataset.zip # removes  zip file to take up less space
        """

rule build_database:
    input:
        "bhv_dataset" # uses input in this folder from download_datasets
    output: # files we will get from makeblastdb 
        "reference/bhv_db.nhr",
        "reference/bhv_db.nin",
        "reference/bhv_db.nsq"
    shell:
        "makeblastdb -in bhv_dataset/ncbi_dataset/data/genomic.fna -out reference/bhv_db -title bhv_db -dbtype nucl" # inputs fasta file with all genomes, makes prefix of output reference/bhv_db, title of database, and it is a nucleotide database

rule blast_query: # running the query with longestb contigs
    input:
        query="assemblies/{sample}/longest_contig.fna", # fasta file with lolngest contig for each sample will be the query 
        db=rules.build_database.output # we will use the output from build_database,  os this rule depends on that rule 
    output:
        results="results/{sample}_blast.csv" # our outptu csv 
    shell:
        """blastn -query {input.query} -db reference/bhv_db -max_hsps 1 -outfmt '6 sacc pident length qstart qend sstart send bitscore evalue stitle' -out {output.results}
        # running nucleotide vs cuncleotide blast, with out wuery of longest contig, uses database we built, limits to highest-scoring segment pair with -max_hsps 1, formats it to be tabelomited, with all the columns we want, and saves it to output 
        """

rule top_5:
    input:
        results="results/{sample}_blast.csv" # uses the output csv file with query results from blast_query
    output:
        "results/{sample}_blast.txt" # # text final formatted correctly with top 5 hits
    shell:
        """
        echo "{wildcards.sample}:" > {output} # sets the sample name as header 
        echo -e "sacc\tpident\tlength\tqstart\tqend\tsstart\tsend\tbitscore\tevalue\tstitle" >> {output} # makes the second line each column name, tab deliminted
        head -n 5 {input} >> {output} # takes the first 4 searchs from blast_query
        """

rule pipeline_report: # concates all the written outsumaries into one file 
    input:
        expand("results/{sample}.report_q2.txt", sample=samples), # from question 2 
        expand("results/{sample}_stats.txt", sample=samples), # from question 3 
        expand("results/{sample}_blast.txt", sample=samples) # from question 3
    output:
        "Sebastian_PipelineReport.txt" # final report 
    shell:
        "cat {input} > {output}" # puts all the input inot the output