#This procedure is to calculate the relative prevalences of the lineages and sublineages of interest by regions and estimate their number of cases

1.- Enter to GISAID (gisaid.org) and download the metadata.
	1.1.- Login
	1.2.- Click on the option "downloads"
	1.3.- In the section "Download packages" download, metadata
	1.4.- In the section "Alignment" download, MSA unmasked

2.- A File named "metadata_tsv_2022_01_22.tar.xz" will be obtained

3.- Decompress the two files to obtain the files: "metadata.tsv"

example command line: 
	tar -xf metadata_tsv_2021_05_25.tar.xz

4.- Use the R script "Preprocessing_metadata_v2.R" to generate an ordered metadata of the complete genomes isolated from Humans

Input file:
	metadata.tsv

example command line: 
	Rscript Preprocessing_metadata_v2.R metadata.tsv

output file:
	ordered_metadata.tsv

5.- Use the R script "Lambda_epidemiological_analysis.R" to calculate the relative prevalences of the sublineages and estimate the number of cases 

Input files:
	ordered_metadata.tsv, positivos_covid.csv (download from https://www.datosabiertos.gob.pe/dataset/casos-positivos-por-covid-19-ministerio-de-salud-minsa), Regions.tsv (file containing the grouping of Peruvian cities in regions), sublineage_SubL1.list, sublineage_SubL2.list (these last files are not specified in the command but have to be in the same folder with the exact names (just changing the number of the sublineage), these files were obtained in section 3_2, step 1)  

example command line: 
	Rscript Lambda_epidemiological_analysis.R ordered_metadata.tsv positivos_covid.csv Regions.tsv

output files:
	Relprev_Lambda_sublineages_reg.jpg, Relprev_Lambda_tot.jpg, Cases_Lambda_sublineages_reg.jpg, Cases_Lambda_tot.jpg, subLambda_epid_data.tsv, subLambda_geno_data.tsv

6.- Use the R script "Gamma_epidemiological_analysis.R" to calculate the relative prevalences of the sublineages and estimate the number of cases 

Input files:
	ordered_metadata.tsv, positivos_covid.csv (download from https://www.datosabiertos.gob.pe/dataset/casos-positivos-por-covid-19-ministerio-de-salud-minsa), Regions.tsv (file containing the grouping of Peruvian cities in regions), sublineage_SubG1.list, sublineage_SubG2.list, sublineage_SubG3.list (these last files are not specified in the command but have to be in the same folder with the exact names (just changing the number of the sublineage), these files were obtained in section 3_2, step 2)

example command line: 
	Rscript Gamma_epidemiological_analysis.R ordered_metadata.tsv positivos_covid.csv Regions.tsv

output files:
	Relprev_Gamma_sublineages_reg.jpg, Relprev_Gamma_tot.jpg, Cases_Gamma_sublineages_reg.jpg, Cases_Gamma_tot.jpg, subGamma_epid_data.tsv, subGamma_geno_data.tsv