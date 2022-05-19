Procedure to obtain sample of genomes from a specific country and a specific VOC between specific dates taking into consideration the regions of isolation of each sample to perform the analysis of transitions of Lambda and Gamma lineages between regions.
*To perform this example is necessary the files obtained in the Preprocessing folder.

1.- Use the files "aligned_290Ns_align_clear_ids.list" and "ordered_metadata.tsv" with the Rscript "sampling_VOC_pg.R" to perform the sampling for the phylogeographic analysis by Peruvian region

example command line:
	To sampling one genome by date by region from Peru of the Omicron variant between 2022-01-01 and 2022-03-03, sampling maximum size is 200:
		Rscript sampling_VOC.R "ordered_metadata.tsv" "aligned_290Ns_align_clear_ids.list" Peru "2022-01-01" "2022-03-03" "s1" "Omicron" 200
	
	To take a random sample of size 200 from Peruvian genomes of the Omicron variant between 2022-01-01 and 2022-03-03, total percentages of genomes by region will be calculated and these proportion wil be considered when taking the sample:		
		Rscript sampling_VOC.R "ordered_metadata.tsv" "aligned_290Ns_align_clear_ids.list" Peru "2022-01-01" "2022-03-03" "sN" "Omicron" 200

	In general:
		Rscript sampling_VOC_pg.R [metadata_file] [seqs_list_file] [country] [date_lower_limit] [date_upper_limit] [type_of_sampling] [VOC] [sample_size]
		
output files:
	Omicron_s1_pg.list
	Omicron_s1_pg.tsv
	Omicron_s1_region_pg.tsv

	or

	Omicron_sN_pg.list
	Omicron_sN_pg.tsv
	Omicron_sN_region_pg.tsv

The file "*.list" contains the list of the IDs of the selected genomes, the "*.tsv" contains the table of the metadata of the selected sequences, and the "*region_pg.tsv" contains the table with the region of each sequence.

2.- Use the script retrieve_sequences_from_txt_notaligned_v2.py to extract the unaligned sequences from the fasta "aligned_290Ns_align_clear.fasta"

example command line:
	python3 retrieve_sequences_from_txt_notaligned_v2.py Omicron_sN.list aligned_290Ns_align_clear.fasta

output file:
	Omicron_sN_aligned.fasta 
	Omicron_sN_seqs.fasta

3.- The file "Omicron_sN_aligned.fasta" contains the sampled sequences aligned, "Omicron_sN_seqs.fasta" contains the sampled sequences not aligned, "Omicron_sN_pg.tsv" (step 1) contains the metadata of the selected sequences.

OPTIONAL:

* The file "Omicron_sN_seqs.fasta" can be uploaded to the web page: https://clades.nextstrain.org/ to remove those that are classified as mediocre or bad by nextstrain:

1.- Wait for the results and click on the button "Download results" in the top-right
2.- Download the "nextclade.csv" file
3.- Use the Rscript "nexclean_filter.R" to create the list of names filtered by nextclade:

example command line:
	Rscript nexclean_filter.R nextclade.csv

output file:
	nexclean.list

4.- Use the script retrieve_sequences_from_txt_notaligned_v2.py to extract the aligned sequences from the previous aligned fasta

example command line: 
	python3 retrieve_sequences_from_txt_notaligned_v2.py nexclean.list Omicron_sN_aligned.fasta

output file:
	nexclean_aligned.fasta nexclean_seqs.fasta


**The files "Omicron_sN_aligned.fasta" (or "nexclean_aligned.fasta") can be used in Beauti to creat the XML files for the BEAST runs.