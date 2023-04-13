Procedure to determine Lambda and Gamma sublineages from a maximum likelihood tree
*this procedure begin with a maximum likelihood tree inferred with iqtree2 (for examples those inferred from the alignments in section 3_1, last step: Gamma_good_seqs.fasta.treefile, Lambda_good_seqs.fasta.treefile)

####Lambda####

1.- Use the script Sublineages_selection_Lambda.R to identify sublineages and take the list of the genomes that belong to each of these sublineages

Input file:
	Lambda_good_seqs.fasta.treefile, nexclean_Peru_Var.tsv (from section 3_1, step 16), Regions.tsv (a file containing the grouping information of the cities by regions) 

example command line:
	Rscript Sublineages_selection_Lambda.R Lambda_good_seqs.fasta.treefile nexclean_Peru_Var.tsv Regions.tsv
	
output file:
	sublineage_SubL1.list, sublineage_SubL2.list, Maximum_likelihood_Lambda.jpg

####Gamma####

2.- Use the script Sublineages_selection_Gamma.R to identify sublineages and take the list of the genomes that belong to each of these sublineages

Input file:
	Gamma_good_seqs.fasta.treefile, nexclean_Peru_Var.tsv (from section 3_1, step 16), Regions.tsv (a file containing the grouping information of the cities by regions) 

example command line:
	Rscript Sublineages_selection_Gamma.R Gamma_good_seqs.fasta.treefile nexclean_Peru_Var.tsv Regions.tsv
	
output file:
	sublineage_SubG1.list, sublineage_SubG2.list, sublineage_SubG3.list, Maximum_likelihood_Gamma.jpg