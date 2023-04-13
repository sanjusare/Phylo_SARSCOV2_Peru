#!usr/bin/env python
import sys
from Bio import SeqIO
from Bio.Seq import Seq
import csv

my_fasta=sys.argv[1]
my_pos=sys.argv[2]
my_out=".".join(my_fasta.split(".")[:-1])+"_spec_sites.fasta"


pos_txt = []
for c in csv.reader(open(my_pos,"r")):
	if c != ['sites']:
            pos_txt.append(int(c[0])-1)
print(pos_txt)

new_aln=[]
for rec in SeqIO.parse(my_fasta,"fasta"):
    filt_seq=[]
    for c in pos_txt:
        filt_seq.extend(rec.seq[c].upper())
    rec.seq=Seq("".join(filt_seq))
    print(rec.seq)
    rec.id=rec.description.split("|")[0]
    rec.description=""	
    new_aln.append(rec)

SeqIO.write(new_aln,my_out,"fasta")
print(my_out)

