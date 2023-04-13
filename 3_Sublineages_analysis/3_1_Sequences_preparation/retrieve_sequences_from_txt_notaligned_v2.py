#!usr/bin/env python

import sys
import csv
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

my_seqs = sys.argv[1]#.txt
my_fasta = sys.argv[2]#.fasta
my_out = my_seqs.split(".")[0]+"_aligned.fasta"
my_out2 = my_seqs.split(".")[0]+"_seqs.fasta"

dic_rec = {}
for rec in SeqIO.parse(my_fasta,"fasta"):
	dic_rec[rec.description] = rec
recs_aligned = []
recs_seqs = []
for c in csv.reader(open(my_seqs,"r")):
    print(c)
    if c[0] in dic_rec:
        recs_aligned.append(dic_rec[c[0]])
        rec_2 = SeqRecord(Seq(str(dic_rec[c[0]].seq).replace("-","")),id=dic_rec[c[0]].id,description=dic_rec[c[0]].description)
        recs_seqs.append(rec_2)
        print(c[0])
SeqIO.write(recs_aligned,my_out,"fasta")
SeqIO.write(recs_seqs,my_out2,"fasta")

			