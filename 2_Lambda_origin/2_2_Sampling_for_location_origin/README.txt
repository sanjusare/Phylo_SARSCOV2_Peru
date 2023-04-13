Procedure to generate the dataset of Lambda genomes to be used in the analysis of the origin of Lambda

1.- Enter to GISAID (gisaid.org) and download the metadata.
	1.1.- Login
	1.2.- Click on the option "downloads"
	1.3.- In the section "Download packages" download, metadata
	1.4.- In the section "Alignment" download, MSA unmasked

2.- A File named "metadata_tsv_2022_01_22.tar.xz" and "sequences_fasta_2021_05_25.tar.xz" will be obtained

3.- Decompress the two files to obtain the files: "metadata.tsv" and "sequences.fasta"

example command line: 
	tar -xf metadata_tsv_2021_05_25.tar.xz
	tar -xf sequences_fasta_2021_05_25.tar.xz

4.- Use the R script "Preprocessing_metadata_v2.R" to generate an ordered metadata of the complete genomes isolated from Humans

Input file:
	metadata.tsv

example command line: 
	Rscript Preprocessing_metadata_v2.R metadata.tsv

output file:
	ordered_metadata.tsv

5.- Use the R script "Sequences_preparation.R" to generate a table with data of the Lambda genomes from selected countries, a list of this genomes and a table to rename the sequences. 

Input files:
	ordered_metadata.tsv

example command line: 
	Rscript Sequences_preparation.R ordered_metadata.tsv

output files:
	sample_initial.tsv, sample_initial.list, sample_rename_initial.tsv

6.- Use the script "retrieve_sequences_from_txt_v3.py" to take the sequences listed in the previously obtained file "sample_initial.list"

Input files:
	sequences.fasta (from step 3)

example command line: 
	python3 retrieve_sequences_from_txt_v3.py sequences.fasta sample_initial.list

output files:
	sample_initial_seqs.fasta

7.- Use the script "renamed_beast_v2.py" to rename the extracted sequences to a more manageable ids

Input files:
	sample_initial_seqs.fasta, sample_rename_initial.tsv (obtained in step 5)

example command line: 
	python3 renamed_beast_v2.py sample_initial_seqs.fasta sample_rename_initial.tsv

output files:
	sample_initial_seqs_renamed.fasta

8.- Align the renamed sequences using ViralMSA.py (https://github.com/niemasd/ViralMSA) to the sequence EPI_ISL_406801 from nucleotide 203 to 29674.

Input files:
	sample_initial_seqs_renamed.fasta

example command line: 
	ViralMSA.py -s sample_initial_seqs_renamed.fasta -r EPI_ISL_406801.fa -e email@gmail.com -o alignment_init -t 5 --omit_ref

output files:
	alignment_init/sample_initial_seqs_renamed.fasta.aln

9.- Remove sequences with more than 290 Ns using the script remove_Ns.py

Input files:
	sample_initial_seqs_renamed.fasta.aln

example command line: 
	python3 remove_Ns.py sample_initial_seqs_renamed.fasta.aln 290

output files:
	sample_initial_seqs_renamed_Ns.fasta

10.- Remove sequences with more than 2 % gaps using the script clear_short_align.py

Input files:
	sample_initial_seqs_renamed_Ns.fasta

example command line: 
	python3 clear_short_align.py sample_initial_seqs_renamed_Ns.fasta 0.02

output files:
	sample_initial_seqs_renamed_Ns_align_clear.fasta

11.- Remove sequences with duplicated IDs using the script remove_duplicate_id.py

Input files:
	sample_initial_seqs_renamed_Ns_align_clear.fasta

example command line: 
	python3 remove_duplicate_id.py sample_initial_seqs_renamed_Ns_align_clear.fasta

output files:
	sample_initial_seqs_renamed_Ns_align_clear_nodup.fasta

12.- Extract IDs from the fasta previously obtained using the script take_fasta_ids.py

Input files:
	sample_initial_seqs_renamed_Ns_align_clear_nodup.fasta

example command line: 
	python3 take_fasta_ids.py sample_initial_seqs_renamed_Ns_align_clear_nodup.fasta

output files:
	sample_initial_seqs_renamed_Ns_align_clear_nodup_ids.list

13.- Extract aligned and not aligned sequences from the fasta obtained in step 11 using the script python3 retrieve_sequences_from_txt_notaligned_v2.py

Input files:
	sample_initial_seqs_renamed_Ns_align_clear.fasta, sample_initial_seqs_renamed_Ns_align_clear_nodup_ids.list

example command line: 
	python3 retrieve_sequences_from_txt_notaligned_v2.py sample_initial_seqs_renamed_Ns_align_clear_nodup_ids.list sample_initial_seqs_renamed_Ns_align_clear_nodup.fasta

output files:
	sample_initial_seqs_renamed_Ns_align_clear_nodup_ids_aligned.fasta, sample_initial_seqs_renamed_Ns_align_clear_nodup_ids_seqs.fasta

14.- Upload the file "sample_initial_seqs_renamed_Ns_align_clear_nodup_ids_seqs.fasta" to the web page: https://clades.nextstrain.org/

15.- Wait for the results and click on the button "Download results" in the top-right

16.- Download the "nextclade.csv" and change its name to "initial.csv"
*optionally you can install nextclade locally to perform the filtering steps (14,15,16) faster.

17.- Use the Rscript "nexclean_lin_v2.R" to create the list of names filtered by nextclade, it is also needed the file sample_initial.tsv from step 5

Input files:
	initial.csv, sample_initial.tsv

example command line: 
	Rscript nexclean_lin_v2.R initial.csv sample_initial.tsv

output files:
	final_metadata_initial.tsv, nexclean_initial.list, location_initial.tsv

18.- Use the Rscript "Lambda_sampling_by_cases.R" to sample Lambda genomes from the filtered genomes according to a range of dates and the approximate number of total genomes you want

Input files:
	final_metadata_initial.tsv (from step 17), Lambda_epidemiological_data.tsv (from section 2_1, step 5)

example command line: 
	Rscript Lambda_sampling_by_cases.R final_metadata_initial.tsv Lambda_epidemiological_data.tsv 9876345 "16" "2020-12-01" "2021-09-30" 200
*third argument informs a seed, fourth informs the name of the sample (arbitrary), fifth and sixth establish the minimum and maximum date to consider, respectively, seventh is the aprroximate number of total genomes to take in the sample.  

output files:
	sample_v2_16.jpg, sample_v2_16.list, sample_v2_16.tsv, Lambda_correlation_plot.jpg
*The first is a figure of the total filtered Lambda genomes (available genomes with high-quality) in grey bars, the number of genomes in the sample in colored bars, the number of genomes that have to be taken to maintain a perfect correlation between genomes and cases in border red bars and the estimated number of Lambda cases in lines. The second is the list of ids of the selected genomes. The third is the metadata of the selected genomes. The fourth is a plot of correlation between Lambda genomes and Lambda cases before and after sampling

19.- Extract the selected sequences from the initial set of genomes (sample_initial_seqs_renamed_Ns_align_clear_nodup_ids_aligned.fasta) obtained in step 13 using the script retrieve_sequences_from_txt_v3.py 

input files:
	sample_initial_seqs_renamed_Ns_align_clear_nodup_ids_aligned.fasta (from step 13), Lsample_v2_16.list

example command line: 
	python3 retrieve_sequences_from_txt_v3.py sample_initial_seqs_renamed_Ns_align_clear_nodup_ids_aligned.fasta sample_v2_16.list

output files:
	sample_v2_16_seqs.fasta

*likelihood mapping can be done in iqtree2 using the following command:
	iqtree2 -s sample_v2_16_seqs.fasta -lmap 10000 -n 0 -m GTR+F+I
*a maximum likelihood tree to obtain the temporal signal can be obtained also in iqtree2 (http://www.iqtree.org/), followed by an analysis in Tempest (http://tree.bio.ed.ac.uk/software/tempest/):
	iqtree2 -s sample_v2_16_seqs.fasta -m GTR+F+I -alrt 1000 -T auto -keep-ident
*this alignment can be also used for Beast2 analysis (http://www.beast2.org/). 