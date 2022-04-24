Procedure to obtain sample of Lambda genomes between specific dates and related genomes of other lineages.
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

2.- Use the script retrieve_sequences_from_txt_notaligned_v2.py to extract the aligned and unaligned sequences from the fasta "aligned_290Ns_align_clear.fasta"

example command line:
	python3 retrieve_sequences_from_txt_notaligned_v2.py Lambda_origin.list aligned_290Ns_align_clear.fasta

output file:
	Lambda_origin_aligned.fasta 
	Lambda_origin_seqs.fasta

3.- Generate a Maximum likelihood tree using Fasttree

example command line:
	FastTree -nt -gtr -gamma -sprlength 1000 -spr 10 -refresh 0.8 -topm 1.5 close 0.75 Lambda_origin_aligned.fasta > Lambda_origin.treefile

output file:
	Lambda_origin.treefile


4.- Use the files "Lambda_origin.tsv" and "Lambda_origin.treefile" with the Rscript "sampling_from_tree_origin.R" to perform the sampling for the Lambda origin estimations

example command line:
		Rscript sampling_from_tree_origin.R "Lambda_origin.tsv" "Lambda_origin.treefile"

	In general:
		Rscript sampling_from_tree_origin.R [metadata_file] [tree_file]
	
output files:
	Lambda_origin_treesample.list
	Lambda_origin_treesample.tsv

The file "*.list" contains the list of the IDs of the selected genomes and the "*.tsv" contains the table of the metadata of the selected sequences.

5.- Use the script retrieve_sequences_from_txt_notaligned_v2.py to extract the unaligned sequences from the fasta "Lambda_origin_aligned.fasta"

example command line:
	python3 retrieve_sequences_from_txt_notaligned_v2.py Lambda_origin_treesample.list Lambda_origin_aligned.fasta

output file:
	Lambda_origin_treesample_aligned.fasta 
	Lambda_origin_treesample_seqs.fasta

6.- The file "Lambda_origin_treesample_aligned.fasta" contains the sampled sequences aligned, "Lambda_origin_treesample_seqs.fasta" contains the sampled sequences not aligned, "Lambda_origin_treesample.tsv" contains the metadata of the selected sequences.

OPTIONAL:

* The file "Lambda_origin_treesample_seqs.fasta" can be uploaded to the web page: https://clades.nextstrain.org/ to remove those that are classified as mediocre or bad by nextstrain:

1.- Wait for the results and click on the button "Download results" in the top-right
2.- Download the "nextclade.csv" file
3.- Use the Rscript "nexclean_filter.R" to create the list of names filtered by nextclade:

example command line:
	Rscript nexclean_filter.R nextclade.csv

output file:
	nexclean.list

4.- Use the script retrieve_sequences_from_txt_notaligned_v2.py to extract the aligned sequences from the previous aligned fasta

example command line: 
	python3 retrieve_sequences_from_txt_notaligned_v2.py nexclean.list Lambda_origin_treesample_aligned.fasta

output file:
	nexclean_aligned.fasta nexclean_seqs.fasta


**The files "Lambda_origin_treesample_aligned.fasta" (or "nexclean_aligned.fasta") can be used in Beauti to creat the XML files for the BEAST runs.