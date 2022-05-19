Procedure to obtain sample of VOC genomes to determine the introductions to Peru
*To perform this example is necessary the files obtained in the Preprocessing folder.

1.- Use the files "aligned_290Ns_align_clear_ids.list" and "ordered_metadata.tsv" with the Rscript "sampling_introduction.R" to perform the sampling for the estimation of the introduction of VOCs to a country

example command line:
	To sampling all the genomes from the VOC of interest (in this example = Gamma) between 45 days before and 15 days after of the first reported genome of the VOC of interest in Peru:
		Rscript sampling_introduction.R "ordered_metadata.tsv" "aligned_290Ns_align_clear_ids.list" "Peru" "Gamma" 45 15
	
	In general:
		Rscript sampling_introductions.R [metadata_file] [seqs_list_file] [country] [variant] [days_before] [days_after]
	
output files:
	Gamma_introduction.list
	Gamma_introduction.tsv

The file "*.list" contains the list of the IDs of the selected genomes and the "*.tsv" contains the table of the metadata of the selected sequences.

2.- Use the script retrieve_sequences_from_txt_notaligned_v2.py to extract the aligned and unaligned sequences from the fasta "aligned_290Ns_align_clear.fasta"

example command line:
	python3 retrieve_sequences_from_txt_notaligned_v2.py Gamma_introduction.list aligned_290Ns_align_clear.fasta

output file:
	Gamma_introduction_aligned.fasta 
	Gamma_introduction_seqs.fasta

3.- Generate a Maximum likelihood tree using Fasttree

example command line:
	FastTree -nt -gtr -gamma -sprlength 1000 -spr 10 -refresh 0.8 -topm 1.5 close 0.75 Gamma_introduction_aligned.fasta > Gamma_introduction.treefile

output file:
	Gamma_introduction.treefile

4.- Use the files "Gamma_introduction.tsv" and "Gamma_introduction.treefile" with the Rscript "sampling_from_tree_introduction.R" to perform the sampling for estimations of the introductions of variants to Peru

example command line:
	Rscript sampling_from_tree_introduction.R "Gamma_introduction.tsv" "Gamma_introduction.treefile"

	In general:
		Rscript sampling_from_tree_introduction.R [metadata_file] [tree_file]
	
output files:
	Gamma_introduction_treesample.list
	Gamma_introduction_treesample.tsv

The file "*.list" contains the list of the IDs of the selected genomes and the "*.tsv" contains the table of the metadata of the selected sequences.

5.- Use the script retrieve_sequences_from_txt_notaligned_v2.py to extract the unaligned sequences from the fasta "Gamma_introduction_aligned.fasta"

example command line:
	python3 retrieve_sequences_from_txt_notaligned_v2.py Gamma_introduction_treesample.list Gamma_introduction_aligned.fasta

output file:
	Gamma_introduction_treesample_aligned.fasta 
	Gamma_introduction_treesample_seqs.fasta

6.- The file "Gamma_introduction_treesample_aligned.fasta" contains the sampled sequences aligned, "Gamma_introduction_treesample_seqs.fasta" contains the sampled sequences not aligned, "Gamma_introduction_treesample.tsv" contains the metadata of the selected sequences.

OPTIONAL:

* The file "Gamma_introduction_treesample_seqs.fasta" can be uploaded to the web page: https://clades.nextstrain.org/ to remove those that are classified as mediocre or bad by nextstrain:

1.- Wait for the results and click on the button "Download results" in the top-right
2.- Download the "nextclade.csv" file
3.- Use the Rscript "nexclean_filter.R" to create the list of names filtered by nextclade:

example command line:
	Rscript nexclean_filter.R nextclade.csv

output file:
	nexclean.list

4.- Use the script retrieve_sequences_from_txt_notaligned_v2.py to extract the aligned sequences from the previous aligned fasta

example command line: 
	python3 retrieve_sequences_from_txt_notaligned_v2.py nexclean.list Gamma_introduction_treesample_aligned.fasta

output file:
	nexclean_aligned.fasta nexclean_seqs.fasta


**The files "Gamma_introduction_treesample_aligned.fasta" (or "nexclean_aligned.fasta") can be used in Beauti to creat the XML files for the BEAST runs.