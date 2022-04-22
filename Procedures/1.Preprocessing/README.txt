Procedure to obtain filtered aligned sequences of SARS-CoV-2 genomes from GISAID (less than 290Ns and less than 2 % gaps in the region corresponding to the nucleotide 203 to 29,674 of the reference genome EPI_ISL_406801).

1.- Enter to GISAID (gisaid.org) and download the sequences and the metadata.
	1.1.- Login
	1.2.- Click on the option "downloads"
	1.3.- In the section "Download packages" download, FASTA and metadata

2.- A File named "sequences_fasta_2022_01_22.tar.xz" and  "metadata_tsv_2022_01_22.tar.xz" will be obtained

3.- Decompress the file to obtain the the files: "metadata.tsv" and "sequences.fasta"

example command line: 
	tar -xf metadata_tsv_2021_05_25.tar.xz

4.- Use the R script "Preprocessing_metadata.R" to generate a table to rename sequences in the fasta, an ordered metadata of the complete genomes isolated from Humans, and a list of the IDs of these genomes

example command line: 
	Rscript Preprocessing_metadata.R metadata.tsv

output files:
	rename_table.tsv
	ordered_metadata.tsv
	selected_seqs.list

5.- Use the python script "retrieve_sequences_from_txt_v3.py" to extract the sequence of the complete genomes isolated from Humans from the fasta decompressed in step 3.

example command line: 
	python3 retrieve_sequences_from_txt_v3.py sequences.fasta selected_seqs.list

output file:
	sequences_seqs.fasta

6.- Use the python script "renamed_beast_v2.py" to rename the sequences in the previous fasta with IDs in a format from where BEAST can easily extract the information.

example command line: 
	python3 renamed_beast_v2.py sequences_seqs.fasta rename_table.tsv

output file:
	sequences_seqs_renamed.fasta

7.- Align the sequences using ViralMSA.py using genome "EPI_ISL_406801.fa" as the reference:

example command line: 
	ViralMSA.py -s sequences_seqs_renamed.fasta -r EPI_ISL_406801.fa -e my_email@gmail -o alignment -t 20 --omit_ref

output file:
	- A folder named alignment will be created with the file named:
		sequences_clear.fasta.aln
	- Change the name to "aligned.fasta"

8.- Remove sequences with more than 290 Ns using the script "remove_Ns.py"

example command line:
	python3 remove_Ns.py aligned.fasta 290

output file:
	aligned_290Ns.fasta

9.- Remove sequences with more than 2 % gaps using the script "clear_short_align.py"

example command line:
	python3 clear_short_align.py VOCs_Peru_aligned_290Ns.fasta 0.02

output file:
	aligned_290Ns_align_clear.fasta

10.- Use the script "take_fasta_ids.py" to take the ids of the last fasta

example command line:
	python3 take_fasta_ids.py VOCs_Peru_aligned_290Ns_align_clear.fasta

output file:
	aligned_290Ns_align_clear_ids.list

11.- Files "aligned_290Ns_align_clear.fasta" and "aligned_290Ns_align_clear_ids.list" contain the filtered sequences in fasta format and a list of their IDs, respectively. Additionally, the file "ordered_metadata.tsv" contains the relevant metadata of all the sequences. These three files are the base for other analysis.