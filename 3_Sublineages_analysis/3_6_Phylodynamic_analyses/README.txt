Procedure to generate figures of the phylodynamics of Lambda and Gamma sublineages based on bayesian analysis

*This procedure begin with combined log files (SubL_BI_fin_v21_seqsf2.log,SubG1_BI_fin_v21_seqsf2.log,SubG2_BI_fin_v21_seqsf2.log,SubG3_BI_fin_v21_seqsf2.log) and tree files with location information (SubL_BI_fin_v21_seqstree_locationf2.trees,SubG1_BI_fin_v21_seqstree_locationf2.trees,SubG2_BI_fin_v21_seqstree_locationf2.trees,SubG3_BI_fin_v21_seqstree_locationf2.trees) of 4 runs of the .xml files created in section 3_5, last step. This procedure also need a maximum clade credibility tree generated from the tree files with location information (SubL_BI_fin_v21_seqstree_location_MCT.trees,SubG1_BI_fin_v21_seqstree_location_MCT.trees,SubG2_BI_fin_v21_seqstree_location_MCT.trees,SubG3_BI_fin_v21_seqstree_location_MCT.trees). The first two files comes from an BEAST2 (http://www.beast2.org/) analysis and the third file is generated using the files SubL_BI_fin_v21_seqstree_locationf2.trees,SubG1_BI_fin_v21_seqstree_locationf2.trees,SubG2_BI_fin_v21_seqstree_locationf2.trees,SubG3_BI_fin_v21_seqstree_locationf2.trees in treeannotator software (https://www.beast2.org/treeannotator/)

1.- Use the R script "lins_origin.R" to generate the figures of date and location of sublineages origin.

Input file:
	SubL_BI_fin_v21_seqstree_location_MCT.trees,SubG1_BI_fin_v21_seqstree_location_MCT.trees,SubG2_BI_fin_v21_seqstree_location_MCT.trees,SubG3_BI_fin_v21_seqstree_location_MCT.trees (those files are not specified as arguments, they have to be in the same folder than the script)

example command line: 
	Rscript lins.origin.R

output file:
	Date_Location_of_Sublineages_origin.jpg

2.- Use the R script "transitions_lins.R" to generate the matrices of transitions. 

Input files:
	SubL_BI_fin_v21_seqstree_locationf2.trees,SubG1_BI_fin_v21_seqstree_locationf2.trees,SubG2_BI_fin_v21_seqstree_locationf2.trees,SubG3_BI_fin_v21_seqstree_locationf2.trees (those files are not specified as arguments, they have to be in the same folder than the script)	

example command line: 
	Rscript transitions_lins.R

output files:
	16_v2tree_location_MCT_Arg.trees, 16_v2tree_location_MCT_Per.trees.

3.- Use the script "clock_rate_sublin.R" to generate the figure of the clock rate distribution

Input files:
	SubL_BI_fin_v21_seqsf2.log,SubG1_BI_fin_v21_seqsf2.log,SubG2_BI_fin_v21_seqsf2.log,SubG3_BI_fin_v21_seqsf2.log (those files are not specified as arguments, they have to be in the same folder than the script)	

example command line: 
	Rscript clock_rate_plot.R

output files:
	Clock_rate_sublineages.jpg