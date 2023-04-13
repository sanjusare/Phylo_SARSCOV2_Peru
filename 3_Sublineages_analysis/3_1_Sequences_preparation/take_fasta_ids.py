from Bio import AlignIO
import pandas as pd
import sys

my_fasta = sys.argv[1]
my_out = ".".join(my_fasta.split(".")[:-1])+"_ids.list"
aln=AlignIO.read(my_fasta,"fasta")


ids=[]
for r in aln:
    ids.append(r.description)

ids_data=pd.DataFrame({"ids":ids})
ids_data.to_csv(my_out,index=False,header=False)

