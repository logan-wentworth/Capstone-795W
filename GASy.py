#! /usr/bin/env python3

#GASy is a comprehensive genome pipeline, that takes a Fastq input and runs it through all
#the steps to push out a polished Fasta file. Some of the parameters are set to a specified
#default to streamline and automate use, while some arguments are accepted for parameters
#that are known to be highly variable. There are options to skip certain portions of the
#program, as it can be used to determine what reference genome to use and may be run again
#using said reference genome. Make sure environment is activated

import subprocess
import urllib.request
import argparse
import os

#Add an option TO run the assmembly portion and removse -s
#Add an option to pull from 16S blast or take manually entered ID and then efetch reference
parser = argparse.ArgumentParser(description="Assemble your genome! Best used with nohup and '&' as this program is time consuming")
parser.add_argument("forward", help="Your forward fastq file to be processed")
parser.add_argument("reverse", help="Your reverse fastq file to be processed")
parser.add_argument("--coverage", "-c", help="Minimum coverage for genome assembly filtering", nargs="?",const="10")
parser.add_argument("--length", "-l", help="Length cutoff for genome assemblyfiltering", nargs="?",const="500")
parser.add_argument("--output", "-o", action="store_true", help="Make a new directory where relevant files will be stored")
parser.add_argument("--plot", "-p", action="store_true", help="Uses blobtool and generates a chart to visualize remaining contigs")
parser.add_argument("--quality", "-q", action="store_true", help="Uses Quast and generates a quality report from filtered genome")
parser.add_argument("--reference", "-r", help="Reference genome if indexing and mapping is desired")
parser.add_argument("--skip", "-s", action="store_true", help="Skips the assembly and annotation if already completed prior")
args = parser.parse_args()

if args.skip is False:
    subprocess.run("trim_scriptV2.sh {0} {1}".format(args.forward, args.reverse), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    subprocess.run("spades.py -1 trimmed-reads/{0} -2 trimmed-reads/{1} -s trimmed-reads/unpaired-{0} -s trimmed-reads/unpaired-{1} -o spades_assembly -t 24".format(args.forward, args.reverse), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    subprocess.run("prokka spades_assembly/contigs.fasta --outdir prokka_output --cpus 24 --mincontiglen 200", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    subprocess.run("grep -o 'product=.*' prokka_output/PROKKA_*.gff | sed 's/product=//g' | sort | uniq -c | sort -nr > protein_content.txt", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    subprocess.run("extract_sequences '16S ribosomal RNA' prokka_output/PROKKA_*.ffn > 16S_sequence.fasta", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    subprocess.run("makeblastdb -in spades_assembly/contigs.fasta -dbtype nucl -out contigs_db", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    subprocess.run("blastn -query 16S_sequence.fasta -db contigs_db -out 16S_vs_contigs.tsv -outfmt 6", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    subprocess.run("blob_blast.sh spades_assembly/contigs.fasta", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
if args.reference is not None:
    subprocess.run("bwa index {0}".format(args.reference), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    subprocess.run("bwa mem -t 24 {0} trimmed-reads/{1} trimmed-reads/{2} > raw_mapped.sam".format(args.reference, args.forward, args.reverse), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    subprocess.run("samtools view -@ 24 -Sb raw_mapped.sam | samtools sort -@ 24 - sorted_mapped", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    subprocess.run("samtools index sorted_mapped.bam", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    subprocess.run("blobtools create -i spades_assembly/contigs.fasta -b sorted_mapped.bam -t contigs.fasta.vs.nt.cul5.1e5.megablast.out -o blob_out", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    subprocess.run("blobtools view -i blob_out.blobDB.json -r all -o blob_taxonomy", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    subprocess.run("grep -v '##' blob_taxonomy.blob_out.blobDB.table.txt | awk -F'\t' '$2 > {0}' | awk -F'\t' '$5 > {1}' | awk -F'\t' '{2}' > contigs_500len_10cov.txt".format(args.length, args.coverage, "{print $1}"), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    subprocess.run("filter_contigs_by_list.py spades_assembly/contigs.fasta contigs_500len_10cov.txt my_new_filtered.fasta", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

if args.quality is True:
    subprocess.run("quast.py my_new_filtered.fasta -o quast_results", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

if args.plot is True:
    subprocess.run("blobtools plot -i blob_out.blobDB.json -r genus", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

if args.output is True:
    subprocess.run("mkdir GASy_output", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    subprocess.run("mv protein_content.txt my_new_filtered.fasta contigs.fasta.vs.nt.cul5.1e5.megablast.out trimmed-reads GASy_output", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if args.quality is True:
        subprocess.run("mv quast_results GASy_output", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if args.plot is True:
        subprocess.run("mv *png GASy_output", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
