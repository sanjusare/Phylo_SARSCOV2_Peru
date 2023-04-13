Procedure to generate figures of the origin of Lambda based on bayesian analysis

*This procedure begin with a combined log files (sample_v2_16f.log) and tree files with location information (16_tree_location_f2.trees) of 4 runs of the .xml files created in section 2_2, last step. This procedure also need a maximum clade credibility tree generated from the tree files with location information (16_tree_location_MCT2.tree). The first two files comes from an BEAST2 (http://www.beast2.org/) analysis and the third file is generated using the file "16_tree_location_f2.trees" in treeannotator software (https://www.beast2.org/treeannotator/)

1.- Use the R script "Phylodynamic_analyses.R" to generate the figures of date and location of Lambda origin, the distributions of trees by date of transitions in both hypotheses and to generate a files of trees where Peru or Argentina is the root 

Input file:
	16_tree_location_f2.trees, 16_tree_location_MCT2.tree, coords.tsv (file manually created with the coordinates of the countries of interes) 

example command line: 
	Rscript Phylodynamic_analyses.R 16_tree_location_MCT2.tree 16_tree_location_f2.trees coords.tsv

output file:
	Date_Location_of_Lambda_origin.jpg, Transitions_distributions.jpg, Hypotheses_world_maps.jpg, 16_v2tree_location_f2_Arg.trees, 16_v2tree_location_f2_Per.trees,
	Summary_of_transitions.tsv

2.- Use tree annotator to generate the MCC tree from the trees when the root was in Argentina or in Peru. 

Input files:
	16_v2tree_location_f2_Arg.trees, 16_v2tree_location_f2_Per.trees

output files:
	16_v2tree_location_MCT_Arg.trees, 16_v2tree_location_MCT_Per.trees.

3.- Use the script "MCC_trees_by_hypothesis.R" to generate the figure of the trees representing the two hypotheses of Lambda origin.

Input files:
	16_v2tree_location_MCT_Arg.trees, 16_v2tree_location_MCT_Per.trees, sample_v2_16.tsv (obtained in section 2_2, step 18), Summary_of_transitions.tsv (from step 1)

example command line: 
	Rscript MCC_trees_by_hypothesis.R 16_v2tree_location_MCT_Arg.trees 16_v2tree_location_MCT_Per.trees sample_v2_16.tsv Summary_of_transitions.tsv

output files:
	MCC_trees_of_Lambda_origin.jpg

4.- Use the script "clock_rate_plot.R" to generate the figure of the clock rate distribution

Input files:
	sample_v2_16f.log

example command line: 
	Rscript clock_rate_plot.R sample_v2_16f.log

output files:
	Clock_rate_distribution.jpg