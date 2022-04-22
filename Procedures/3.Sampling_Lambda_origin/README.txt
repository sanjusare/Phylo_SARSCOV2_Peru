Procedure to obtain sample of genomes from a specific country and a specific VOC between specific dates.
*To perform this example is necessary the files obtained in the Preprocessing folder.

1.- Use the files "aligned_290Ns_align_clear_ids.list" and "ordered_metadata.tsv" with the Rscript "sampling_Lambda_origin.R" to perform the sampling for the Lambda origin estimation

example command line:
	To sampling all the Lambda genomes, all the Peruvian genomes from other lineages, and one of each day of each lineage of each other country between 2020-06-01 and 2021-31-01:
		Rscript sampling_Lambda_origin.R "ordered_metadata.tsv" "aligned_290Ns_align_clear_ids.list" "2022-06-01" "2021-31-01"
	
	In general:
		Rscript sampling_Lambda_origin.R [metadata_file] [seqs_list_file] [date_lower_limit] [date_upper_limit]
	
output files:
	Lambda_origin.list
	Lambda_origin.tsv

The file "*.list" contains the list of the IDs of the selected genomes and the "*.tsv" contains the table of the metadata of the selected sequences.

2.- Use the script retrieve_sequences_from_txt_notaligned_v2.py to extract the unaligned sequences from the fasta "aligned_290Ns_align_clear.fasta"

example command line:
	python3 retrieve_sequences_from_txt_notaligned_v2.py Lambda_origin.list aligned_290Ns_align_clear.fasta

output file:
	Lambda_origin_aligned.fasta 
	Lambda_origin_seqs.fasta

3.- Generate a Maximum likelihood tree using Fasttree

example command line:
	FastTree -nt -gtr -gamma -sprlength 1000 -spr 10 -refresh 0.8 -topm 1.5 close 0.75 Lambda_origin_aligned.fasta > Lambda_origin_aligned.treefile

output file:
	Lambda_origin_aligned.treefile


