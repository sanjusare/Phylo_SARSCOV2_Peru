Procedure to generate figures of general epidemiological data (Figure 1).

1.- Enter to GISAID (gisaid.org) and download the metadata.
	1.1.- Login
	1.2.- Click on the option "downloads"
	1.3.- In the section "Download packages" download, metadata

2.- A File named "metadata_tsv_2022_01_22.tar.xz" will be obtained

3.- Decompress the file to obtain the the files: "metadata.tsv"

example command line: 
	tar -xf metadata_tsv_2021_05_25.tar.xz

4.- Use the R script "Preprocessing_metadata.R" to generate an ordered metadata of the complete genomes isolated from Humans

Input file:
	metadata.tsv

example command line: 
	Rscript Preprocessing_metadata.R metadata.tsv

output file:
	ordered_metadata.tsv

5.- Use the R script "Epidemiological_Figures.R" to generate Figures of Cases and Deaths, Stringency index and Vaccinations, and Relative Prevalences of VOCIs.

Input files:
	ordered_metadata.tsv, owid-covid-data.csv (#download from https://github.com/owid/covid-19-data/tree/master/public/data)

example command line: 
	Rscript Epidemiological_figures.R ordered_metadata.tsv owid-covid-data.csv positivos_covid.csv

output files:
	Reported_cases_deaths.jpg, Stringency_index_Vaccinations.jpg, Relative_Prevalence_VOCI.jpg, General_analysis.jpg, Reported_deaths_cases_country.jpg, Test_by_inhabitant.jpg, Reported_deaths_cases_tests_country.jpg   


