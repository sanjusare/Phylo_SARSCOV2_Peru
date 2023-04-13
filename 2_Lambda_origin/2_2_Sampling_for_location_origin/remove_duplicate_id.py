from Bio import SeqIO
import sys

my_fasta=sys.argv[1]
my_out = my_fasta.split(".")[0]+"_nodup."+my_fasta.split(".")[1]


with open(my_out, 'a') as outFile:
    record_ids = list()
    for record in SeqIO.parse(my_fasta, 'fasta'):
        if record.id not in record_ids:
            record_ids.append(record.id)
            SeqIO.write(record, outFile, 'fasta')