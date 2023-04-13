#This procedure is to sampling the genomes that will be used for the inferences of the phylodynamics inside Peru of the lineages and sublineages of interest

####Lambda####

1.- Use the R script "subLambda_sampling.R" to obtain the samples of genomes to perform the bayesian analyses and a plot of the correlation between cases and genomes before and after sampling

Input files:
	subLambda_epid_data.tsv, subLambda_geno_data.tsv (both from section 3_4, step 5)

example command line: 
	Rscript subLambda_sampling.R subLambda_epid_data.tsv subLambda_geno_data.tsv

output files:
	Lambda_sampling_correlation.jpg, Lambda_BI_fin_v2.list, Lambda_region_fin_v2.tsv, Lambda_meta_fin_v2.tsv  

2.- Use the script "retrieve_sequences_from_txt_v3.py" to obtain the genomes from the file containing all the Lambda high-quality sequences from Peru (Lambda_good_seqs.fasta)

Input files:
	Lambda_BI_fin_v2.list, Lambda_good_seqs.fasta (section 3_1, step 18)

example command line: 
	python3 retrieve_sequences_from_txt_v3.py Lambda_good_seqs.fasta Lambda_BI_fin_v2.list

output files:
	Lambda_BI_fin_v2_seqs.fasta

*This last file can be used for maximum likelihood phylogenetic inference in iqtree2 (http://www.iqtree.org/) with the following command:
	iqtree2 -s Lambda_BI_fin_v2_seqs.fasta -m GTR+F+I -B 1000 -alrt 1000 -T auto -keep-ident
*likelihood mapping can be done in iqtree2 using the following command:
	iqtree2 -s Lambda_BI_fin_v2_seqs.fasta -lmap 10000 -n 0 -m GTR+F+I
*the temporal signal can be obtained from the maximum likelihood trees using Tempest (http://tree.bio.ed.ac.uk/software/tempest/) 
*The file Lambda_BI_fin_v2_seqs.fasta together with the file Lambda_region_fin_v2.tsv were used to generate the .xml files for bayesian analysis. 

####Gamma#### 

3.- Use the R script "subGamma_sampling.R" to obtain the samples of genomes to perform the bayesian analyses and a plot of the correlation between cases and genomes before and after sampling

Input files:
	subGamma_epid_data.tsv, subGamma_geno_data.tsv (both from section 3_4, step 6)

example command line: 
	Rscript subGamma_sampling.R subGamma_epid_data.tsv subGamma_geno_data.tsv

output files:
	Gamma_sampling_correlation.jpg, SubG1_BI_fin_v2.list, SubG1_region_fin_v2.tsv, SubG1_meta_fin_v2.tsv, SubG2_BI_fin_v2.list, SubG2_region_fin_v2.tsv, SubG2_meta_fin_v2.tsv, SubG3_BI_fin_v2.list, SubG3_region_fin_v2.tsv, SubG3_meta_fin_v2.tsv 

4.- Use the script "retrieve_sequences_from_txt_v3.py" to obtain the genomes from the file containing all the Gamma high-quality sequences from Peru (Gamma_good_seqs.fasta)

Input files:
	SubG1_BI_fin_v2.list, SubG2_BI_fin_v2.list, SubG3_BI_fin_v2.list, Gamma_good_seqs.fasta (section 3_1, step 18)

example command line: 
	python3 retrieve_sequences_from_txt_v3.py Gamma_good_seqs.fasta SubG1_BI_fin_v2.list
	python3 retrieve_sequences_from_txt_v3.py Gamma_good_seqs.fasta SubG2_BI_fin_v2.list
	python3 retrieve_sequences_from_txt_v3.py Gamma_good_seqs.fasta SubG3_BI_fin_v2.list

output files:
	SubG1_BI_fin_v2_seqs.fasta, SubG2_BI_fin_v2_seqs.fasta, SubG3_BI_fin_v2_seqs.fasta

*These last files can be used for maximum likelihood phylogenetic inference in iqtree2 (http://www.iqtree.org/) with the following command:
	iqtree2 -s SubG1_BI_fin_v2_seqs.fasta -m GTR+F+I -B 1000 -alrt 1000 -T auto -keep-ident
	iqtree2 -s SubG2_BI_fin_v2_seqs.fasta -m GTR+F+I -B 1000 -alrt 1000 -T auto -keep-ident
	iqtree2 -s SubG3_BI_fin_v2_seqs.fasta -m GTR+F+I -B 1000 -alrt 1000 -T auto -keep-ident
*likelihood mapping can be done in iqtree2 using the following command:
	iqtree2 -s SubG1_BI_fin_v2_seqs.fasta -lmap 10000 -n 0 -m GTR+F+I
	iqtree2 -s SubG2_BI_fin_v2_seqs.fasta -lmap 10000 -n 0 -m GTR+F+I
	iqtree2 -s SubG3_BI_fin_v2_seqs.fasta -lmap 10000 -n 0 -m GTR+F+I
*the temporal signal can be obtained from the maximum likelihood trees using Tempest (http://tree.bio.ed.ac.uk/software/tempest/) 
*The files SubG1_BI_fin_v2_seqs.fasta, SubG2_BI_fin_v2_seqs.fasta, SubG3_BI_fin_v2_seqs.fasta together with the files SubG1_region_fin_v2.tsv, SubG2_region_fin_v2.tsv, SubG3_region_fin_v2.tsv were used to generate the .xml files for bayesian analysis. 

