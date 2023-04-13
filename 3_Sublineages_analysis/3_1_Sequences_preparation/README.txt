Procedure to prepare the Lambda and Gamma genomes to use in the maximum likelihood analysis

1.- Enter to GISAID (gisaid.org) and download the metadata.
	1.1.- Login
	1.2.- Click on the option "downloads"
	1.3.- In the section "Download packages" download, metadata
	1.4.- Download all the Peruvian genomes from GISAID

2.- A File named "metadata_tsv_2022_01_22.tar.xz" will be obtained, and a file containing the Peruvian fasta sequences that we will named VOCs_Peru.fasta

3.- Decompress the file metadata_tsv_2022_01_22.tar.xz to obtain the file: "metadata.tsv"

example command line: 
	tar -xf metadata_tsv_2021_05_25.tar.xz

4.- Use the R script "Preprocessing_metadata_v3.R" to generate an ordered metadata of the complete genomes isolated from Humans and a table to rename fasta to a more manageable names

Input file:
	metadata.tsv

example command line: 
	Rscript Preprocessing_metadata_v3.R metadata.tsv

output file:
	ordered_metadata.tsv, rename_table.tsv

5.- Use the python script "renamed_beast_v2.py" to rename the VOCs_Peru.fasta

Input file:
	VOCs_Peru.fasta, rename_table.tsv

example command line: 
	python3 renamed_beast_v2.py VOCs_Peru.fasta rename_table.tsv

output file:
	VOCs_Peru_renamed.fasta

6.- Align the renamed sequences using ViralMSA.py (https://github.com/niemasd/ViralMSA) to the sequence EPI_ISL_406801 from nucleotide 203 to 29674.

Input files:
	VOCs_Peru_renamed.fasta

example command line: 
	ViralMSA.py -s VOCs_Peru_renamed.fasta -r EPI_ISL_406801.fa -e email@gmail.com -o alignment_init -t 5 --omit_ref

output files:
	alignment_init/VOCs_Peru_renamed.fasta.aln

7.- Remove sequences with more than 290 Ns using the script remove_Ns.py

Input files:
	VOCs_Peru_renamed.fasta.aln

example command line: 
	python3 remove_Ns.py VOCs_Peru_renamed.fasta.aln 290

output files:
	VOCs_Peru_renamed_Ns.fasta

8.- Remove sequences with more than 2 % gaps using the script clear_short_align.py

Input files:
	VOCs_Peru_renamed_Ns.fasta

example command line: 
	python3 clear_short_align.py VOCs_Peru_renamed_Ns.fasta 0.02

output files:
	VOCs_Peru_renamed_Ns_align_clear.fasta

9.- Use the script "take_fasta_ids.py" to take the ids of the previous fasta

Input files:
	VOCs_Peru_renamed_Ns_align_clear.fasta

example command line: 
	python3 take_fasta_ids.py VOCs_Peru_renamed_Ns_align_clear.fasta

output files:
	VOCs_Peru_renamed_Ns_align_clear_ids.list

10.- Use the script "retrieve_sequences_from_txt_notaligned_v2.py" to extract the unaligned sequences from the previous fasta VOCs_Peru_renamed_Ns_align_clear.fasta (from step 8)

Input files:
	VOCs_Peru_renamed_Ns_align_clear.fasta (step 8), VOCs_Peru_renamed_Ns_align_clear_ids.list

example command line: 
	python3 retrieve_sequences_from_txt_notaligned_v2.py VOCs_Peru_renamed_Ns_align_clear_ids.list VOCs_Peru_renamed_Ns_align_clear.fasta

output files:
	VOCs_Peru_renamed_Ns_align_clear_ids_aligned.fasta, VOCs_Peru_renamed_Ns_align_clear_ids_seqs.fasta

11.- Upload the file "VOCs_Peru_renamed_Ns_align_clear_ids_seqs.fasta" to the web page: https://clades.nextstrain.org/

12.- Wait for the results and click on the button "Download results" in the top-right

13.- Download the "nextclade.csv"
*optionally you can install nextclade locally to perform the filtering steps (11,12,13) faster.

14.- Use the script "nexclean_fRt.R" to create the list of names filtered by nextclade:

Input files:
	nextclade.csv

example command line: 
	Rscript nexclean_fRt.R nextclade.csv

output files:
	nexclean.list

15.- Use the script "retrieve_sequences_from_txt_notaligned_v2.py" to extract the aligned sequences from the previous aligned fasta (from step 10):

Input files:
	VOCs_Peru_renamed_Ns_align_clear_ids_aligned.fasta, nexclean.list

example command line: 
	python3 retrieve_sequences_from_txt_notaligned_v2.py nexclean.list VOCs_Peru_renamed_Ns_align_clear_ids_aligned.fasta

output files:
	nexclean_aligned.fasta, nexclean_seqs.fasta

16.- Use the file "nexclean.list" (step 14) and the Rscript "all_metadata_fRt.R" to extract just the filtered sequences from the metadata (step 4, ordered_metadata):

Input files:
	ordered_metadata.tsv, nexclean.list

example command line: 
	Rscript all_metadata_fRt.R ordered_metadata.tsv nexclean.list

output files:
	nexclean_Peru_Var.tsv

17.- Use the Rscript "variants.R" to extract the a list of Gamma and Lambda genomes from Peru metadata (step 16, nexclean_Peru_Var.tsv):

Input files:
	nexclean_Peru_Var.tsv

example command line: 
	Rscript variants.R nexclean_Peru_Var.tsv

output files:
	Lambda_good.list, Gamma_good.list	

18.- Use the script "retrieve_sequences_from_txt_v3.py" to extract the aligned sequences from the fasta nexclean_aligned.fasta (from step 15)

Input files:
	nexclean_aligned.fasta, Lambda_good.list, Gamma_good.list

example command line: 
	python3 retrieve_sequences_from_txt_v3.py nexclean_aligned.fasta Lambda_good.list
	python3 retrieve_sequences_from_txt_v3.py nexclean_aligned.fasta Gamma_good.list

output files:
	Lambda_good_seqs.fasta, Gamma_good_seqs.fasta
	
*These last files can be used for maximum likelihood phylogenetic inference in iqtree2 (http://www.iqtree.org/) with the following command:
	iqtree2 -s Lambda_good_seqs.fasta -m GTR+F+I -B 1000 -alrt 1000 -T auto -keep-ident
	iqtree2 -s Gamma_good_seqs.fasta -m GTR+F+I -B 1000 -alrt 1000 -T auto -keep-ident
*likelihood mapping can be done in iqtree2 using the following command:
	iqtree2 -s Lambda_good_seqs.fasta -lmap 10000 -n 0 -m GTR+F+I
	iqtree2 -s Gamma_good_seqs.fasta -lmap 10000 -n 0 -m GTR+F+I
*the temporal signal can be obtained from the maximum likelihood trees using Tempest (http://tree.bio.ed.ac.uk/software/tempest/)


