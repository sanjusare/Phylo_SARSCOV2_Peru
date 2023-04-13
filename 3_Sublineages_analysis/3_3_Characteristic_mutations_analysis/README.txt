Procedure to find characteristic mutations that represent each of the identified sublineages
*this procedure begin with the list of the tips of the sublineages created in section 3_2

####Lambda####

1.- Use the script retrieve_sequences_from_txt_v3.py to take the aligned sequences corresponding to each sublineage

Input files:
	Lambda_good_seqs.fasta (from section 3_1, step 18), sublineage_SubL1.list, sublineage_SubL2.list (both from section 3_2, step 1).

example command line:
	python3 retrieve_sequences_from_txt_v3.py Lambda_good_seqs.fasta sublineage_SubL1.list
	python3 retrieve_sequences_from_txt_v3.py Lambda_good_seqs.fasta sublineage_SubL2.list
	
output files:
	sublineage_SubL1_seqs.fasta, sublineage_SubL2_seqs.fasta

2.- Use the script frequencies_v2.py to calculate the relative frequencies of each base in each position of the alignment

Input files:
	sublineage_SubL1_seqs.fasta, sublineage_SubL2_seqs.fasta 

example command line:
	python3 frequencies_v2.py sublineage_SubL1_seqs.fasta
	python3 frequencies_v2.py sublineage_SubL2_seqs.fasta

output files:
	sublineage_SubL1_seqs.fasta.freq, sublineage_SubL2_seqs.fasta.freq

3.- Use the script Mutations_frequency_analysis_Lambda.R to identify mutations with more than 0.8 relative frequency in one lineage and less than 0.2 in other sublineages

Input files:
	sublineage_SubL1_seqs.fasta.freq, sublineage_SubL2_seqs.fasta.freq

example command line:
	Rscript Mutations_frequency_analysis.R sublineage_SubL1_seqs.fasta.freq sublineage_SubL2_seqs.fasta.freq

output files:
	Characteristic_subLambda_mutations.xlsx, selected_Lambda_position.list

4.- Use the script extr_subalign.py to extract the positions of interest from the alignment of the Lambda sequences (from section 3_1, step 18)

Input files:
	Lambda_good_seqs.fasta, selected_Lambda_position.list

example command line:
	python3 extr_subalign.py Lambda_good_seqs.fasta selected_Lambda_position.list

output files:
	Lambda_good_seqs_spec_sites.fasta

5.- Use the script Sublineages_msa_Lambda.R to generate the figure of the tree associated with an alignment of the identified positions

Input files:
	Lambda_good_seqs_spec_sites.fasta, nexclean_Peru_Var.tsv (from section 3_1, step 16), Lambda_good_seqs.fasta.treefile (obtained in section 3_1 last step) 

example command line:
	Rscript Sublineages_msa_Lambda.R Lambda_good_seqs.fasta.treefile nexclean_Peru_Var.tsv Lambda_good_seqs_spec_sites.fasta

output files:
	msa_tree_Lambda.jpg


####Gamma####

6.- Use the script retrieve_sequences_from_txt_v3.py to take the aligned sequences corresponding to each sublineage

Input files:
	Gamma_good_seqs.fasta (from section 3_1, step 18), sublineage_SubG1.list, sublineage_SubG2.list, sublineage_SubG3.list (all lists from section 3_2, step 2).

example command line:
	python3 retrieve_sequences_from_txt_v3.py Gamma_good_seqs.fasta sublineage_SubG1.list
	python3 retrieve_sequences_from_txt_v3.py Gamma_good_seqs.fasta sublineage_SubG2.list
	python3 retrieve_sequences_from_txt_v3.py Gamma_good_seqs.fasta sublineage_SubG3.list
	
output files:
	sublineage_SubG1_seqs.fasta, sublineage_SubG2_seqs.fasta, sublineage_SubG3_seqs.fasta

7.- Use the script frequencies_v2.py to calculate the relative frequencies of each base in each position of the alignment

Input files:
	sublineage_SubG1_seqs.fasta, sublineage_SubG2_seqs.fasta, sublineage_SubG3_seqs.fasta

example command line:
	python3 frequencies_v2.py sublineage_SubG1_seqs.fasta
	python3 frequencies_v2.py sublineage_SubG2_seqs.fasta
	python3 frequencies_v2.py sublineage_SubG3_seqs.fasta

output files:
	sublineage_SubG1_seqs.fasta.freq, sublineage_SubG2_seqs.fasta.freq, sublineage_SubG3_seqs.fasta.freq

8.- Use the script Mutations_frequency_analysis_Gamma.R to identify mutations with more than 0.8 relative frequency in one lineage and least than 0.2 in other sublineages

Input files:
	sublineage_SubG1_seqs.fasta.freq, sublineage_SubG2_seqs.fasta.freq, sublineage_SubG3_seqs.fasta.freq

example command line:
	Rscript Mutations_frequency_analysis_Gamma.R sublineage_SubG1_seqs.fasta.freq sublineage_SubG2_seqs.fasta.freq sublineage_SubG3_seqs.fasta.freq

output files:
	Characteristic_subGamma_mutations.xlsx, selected_Gamma_position.list

9.- Use the script extr_subalign.py to extract the positions of interest from the alignment of the Gamma sequences (from section 3_1, step 18)

Input files:
	Gamma_good_seqs.fasta, selected_Gamma_position.list

example command line:
	python3 extr_subalign.py Gamma_good_seqs.fasta selected_Gamma_position.list

output files:
	amma_good_seqs_spec_sites.fasta

5.- Use the script Sublineages_msa_Gamma.R to generate the figure of the tree associated with an alignment of the identified positions

Input files:
	Gamma_good_seqs_spec_sites.fasta, nexclean_Peru_Var.tsv (from section 3_1, step 16), Gamma_good_seqs.fasta.treefile (obtained in section 3_1 last step) 

example command line:
	Rscript Sublineages_msa_Gamma.R Gamma_good_seqs.fasta.treefile nexclean_Peru_Var.tsv Gamma_good_seqs_spec_sites.fasta

output files:
	msa_tree_Gamma.jpg
