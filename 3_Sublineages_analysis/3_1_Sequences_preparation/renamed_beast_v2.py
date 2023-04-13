#!usr/bin/env python

import sys
import csv
from Bio import SeqIO
#import pandas as pd

my_fasta = sys.argv[1]
my_table = sys.argv[2]
my_out = my_fasta.split(".")[0]+"_renamed."+my_fasta.split(".")[1]

#names=pd.read_csv(my_table,sep="\t",header=None,dtype=str)

# making dic of table
dic_table = {}
for row in csv.reader(open(my_table),delimiter='\t'):
    dic_table[row[0]] = row[1]

count=0    
rec_new_id = []
for r in SeqIO.parse(my_fasta,"fasta"):
    if r.description in dic_table:
        r.id = dic_table[r.description]
        r.description=""
        rec_new_id.append(r)
        count=count+1
        print(count) 
SeqIO.write(rec_new_id,my_out,"fasta")
print(my_out)



'''
rec_new_id = []
for r in SeqIO.parse(my_fasta,"fasta"):
    row=names.loc[names[0]==r.description]
    r.id = row.iat[0,1]
    r.description = ""
    rec_new_id.append(r)

SeqIO.write(rec_new_id,my_out,"fasta")
print(my_out)
'''
