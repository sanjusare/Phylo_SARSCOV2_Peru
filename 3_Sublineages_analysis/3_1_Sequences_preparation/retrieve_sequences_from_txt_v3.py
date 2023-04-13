#!usr/bin/env python

import sys
import csv
from Bio import SeqIO

my_fasta = sys.argv[1]#fasta
my_seqs = sys.argv[2]#list
my_out = my_seqs.split(".")[0]+"_seqs.fasta"

seqs_txt = []
for c in csv.reader(open(my_seqs,"r")):
    print(c)
    seqs_txt.extend(c)
    
def parse_fasta(fname, fseqs):
    rec_seqs=[]
    for rec in SeqIO.parse(fname,"fasta"):
        if rec.description in seqs_txt:
            rec_seqs.append(rec)
            yield rec
    return(rec_seqs)
    
output=parse_fasta(my_fasta, my_seqs)

SeqIO.write(output,my_out,"fasta")

