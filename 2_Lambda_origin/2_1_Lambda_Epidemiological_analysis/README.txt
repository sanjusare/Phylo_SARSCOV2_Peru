Procedure to generate figures of Lambda epidemiological data (Figure 2,S1-S4).

1.- Enter to GISAID (gisaid.org) and download the metadata.
	1.1.- Login
	1.2.- Click on the option "downloads"
	1.3.- In the section "Download packages" download, metadata

2.- A File named "metadata_tsv_2022_01_22.tar.xz" will be obtained

3.- Decompress the file to obtain the the files: "metadata.tsv"

example command line: 
	tar -xf metadata_tsv_2021_05_25.tar.xz

4.- Use the R script "Preprocessing_metadata.R" to generate an ordered metadata (the same obtained in 1_General_Epidemiological_analyses) of the complete genomes isolated from Humans

Input file:
	metadata.tsv

example command line: 
	Rscript Preprocessing_metadata.R metadata.tsv

output file:
	ordered_metadata.tsv

5.- Use the R script "Lambda_epidemiological_Figures.R" to generate Figures of first Lambda genomes, Lambda prevalence and number of genomes by cities, fitted relative prevalences, estimated cases, Lambda genomes and sampling proportions. Also, it will create a table containing the information used in the plots.

Input files:
	ordered_metadata.tsv, owid-covid-data.csv (#download from https://github.com/owid/covid-19-data/tree/master/public/data)

example command line: 
	Rscript Lambda_epidemiological_figures.R ordered_metadata.tsv owid-covid-data.csv

output files:
	First_Lambda_genomes.jpg, Cities_Lambda_relative_prevalence.jpg, Cities_Lambda_genomes.jpg, Fitted_relprev_lam_cases_country.jpg, Lambda_genomes_sampprop_country.jpg, Lambda_epidemiological_data.tsv


