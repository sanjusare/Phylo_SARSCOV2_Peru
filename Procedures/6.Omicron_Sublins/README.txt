Procedure to obtain sample of genomes from a specific country and a specific VOC to perform an analysis of new lineages in that country.
*To perform this example is necessary the files obtained in the Preprocessing folder.

1.- Use the files "aligned_290Ns_align_clear_ids.list" and "ordered_metadata.tsv" with the Rscript "sampling_lineages_id.R" to perform the sampling for the identification of new lineages of an specific VOCI

example command line:
	To sampling all the genomes from Peru of the VOCI selected and one genome of each lineage of each available country of the VOCI analyzed:
		Rscript sampling_lineages_id.R "ordered_metadata.tsv" "aligned_290Ns_align_clear_ids.list" Peru "Omicron"

	In general:
		Rscript sampling_lineages_id.R [metadata_file] [seqs_list_file] [country] [VOC]
		
output files:
	Omicron_Peru_sample_lin_id.list
	Omicron_Peru_sample_lin_id.tsv

The file "*.list" contains the list of the IDs of the selected genomes, the "*.tsv" contains the table of the metadata of the selected sequences.

2.- Use the script retrieve_sequences_from_txt_notaligned_v2.py to extract the unaligned sequences from the fasta "aligned_290Ns_align_clear.fasta"

example command line:
	python3 retrieve_sequences_from_txt_notaligned_v2.py Omicron_Peru_sample_lin_id.list aligned_290Ns_align_clear.fasta

output file:
	Omicron_Peru_sample_lin_id_aligned.fasta Omicron_Peru_sample_lin_id_seqs.fasta

3.- Generate a Maximum likelihood tree using IQtree

example command line:
	iqtree2 -s Omicron_Peru_sample_lin_id_aligned.fasta -m GTR+F+I -alrt 1000 -B 1000 -T auto -keep-ident

output file:
	Omicron_Peru_sample_lin_id_aligned.fasta.treefile


4.- Use the files "Omicron_Peru_sample_lin_id.tsv" and "Omicron_Peru_sample_lin_id_aligned.fasta.treefile" with the Rscript "select_new_lins_from_tree.R" to perform the identification and extractions of the tips that correspond to the identified lineages 

example command line:
	Rscript select_new_lins_from_tree.R "Omicron_Peru_sample_lin_id.tsv" "Omicron_Peru_sample_lin_id_aligned.fasta.treefile" Peru "Omicron"

	In general:
		Rscript select_new_lins_from_tree.R [metadata_file] [tree_file] [Country] [VOC]
	
output files:
	Peru_Omicron_new.list
	Peru_Omicron_new.tsv
	Peru_Omicron1_new.list
	Peru_Omicron1_new.tsv
	Peru_Omicron2_new.list
	Peru_Omicron2_new.tsv

The file "*.list" contains the list of the IDs of the selected genomes and the "*.tsv" contains the table of the metadata of the selected sequences. Peru_Omicron_new contains the information of just the genomes corresponding to the identified lineage, Peru_Omicron_new1 contains those genomes plus the genomes below one node back and Peru_Omicron_new2 contains the genomes of the identified lineages plus the genomes below two nodes back.

5.- Use the script retrieve_sequences_from_txt_notaligned_v2.py to extract the unaligned sequences from the fasta "Omicron_Peru_sample_lin_id_aligned.fasta"

example command line:
	python3 retrieve_sequences_from_txt_notaligned_v2.py Peru_Omicron_new.list Omicron_Peru_sample_lin_id_aligned.fasta

output file:
	Peru_Omicron_new_aligned.fasta 
	Peru_Omicron_new_seqs.fasta

6.- The file "Peru_Omicron_new_aligned.fasta.fasta" contains the sampled sequences aligned, "Peru_Omicron_new_seqs.fasta" contains the sampled sequences not aligned, "Peru_Omicron_new.tsv" contains the metadata of the selected sequences.

OPTIONAL:

* The file "Peru_Omicron_new_seqs.fasta" can be uploaded to the web page: https://clades.nextstrain.org/ to remove those that are classified as mediocre or bad by nextstrain:

1.- Wait for the results and click on the button "Download results" in the top-right
2.- Download the "nextclade.csv" file
3.- Use the Rscript "nexclean_filter.R" to create the list of names filtered by nextclade:

example command line:
	Rscript nexclean_filter.R nextclade.csv

output file:
	nexclean.list

4.- Use the script retrieve_sequences_from_txt_notaligned_v2.py to extract the aligned sequences from the previous aligned fasta

example command line: 
	python3 retrieve_sequences_from_txt_notaligned_v2.py nexclean.list Peru_Omicron_new_aligned.fasta

output file:
	nexclean_aligned.fasta nexclean_seqs.fasta


**The files "Peru_Omicron_new_aligned.fasta" (or "nexclean_aligned.fasta") can be used in Beauti to create the XML files for the BEAST runs.