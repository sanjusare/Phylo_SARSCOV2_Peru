# 1:set working directory####
args <- commandArgs(trailingOnly = TRUE)
# 2:load libraries and establish settings####
library(treeio)
library(ggtree)
library(dplyr)
library(ggplot2)
library(ggnewscale)
library(ggpubr)
theme_set(theme_bw(base_size = 12))
Sys.setlocale("LC_TIME", "English")
collin <- c(SubG1 = "dodgerblue", SubG2 = "firebrick", SubG3 = "#ad71ab")
# 3:read files####
args <- c("Gamma_good_seqs.fasta.treefile", "nexclean_Peru_Var.tsv",
    "Gamma_good_seqs_spec_sites.fasta")
# read the maximum likelihoood tree
t <- read.iqtree(args[1])
# read the metadata of genomes from Peru and filter those
# that are present in the tree
d <- read.csv(args[2], sep = "\t") %>%
    dplyr::filter(name_analysis %in% t@phylo[["tip.label"]]) %>%
    mutate(date = as.Date(date))
# read the fasta of the selected positions
sequences <- read.fasta(args[3])
# 4:Sublineages identification#### selection of sublineages
p_t <- ggtree(t, color = "lightgrey") %<+% d
r_node <- p_t[["data"]]$parent[p_t[["data"]]$node == p_t[["data"]]$parent]
nodes <- p_t[["data"]] %>%
    dplyr::filter(isTip == FALSE & node != r_node & SH_aLRT >=
        70 & UFboot >= 70)
sel_nodes <- c()
tips <- list()
for (i in 1:nrow(nodes)) {
    print(i)
    t1 <- tree_subset(t, node = nodes$node[i], levels_back = 0)
    if (length(t1@phylo[["tip.label"]]) >= 400) {
        tips[[as.character(nodes$node[i])]] <- t1@phylo[["tip.label"]]
        sel_nodes <- append(sel_nodes, nodes$node[i])
    }
}
sel_lins <- list()
for (i in 1:length(sel_nodes)) {
    print(i)
    if (i != length(sel_nodes)) {
        test_tips <- tips[[i]][!(tips[[i]] %in% tips[[i + 1]])]
    } else {
        test_tips <- tips[[i]]
    }
    if (length(test_tips) >= 400) {
        sel_lins[[as.character(sel_nodes[i])]] <- test_tips
    }
}
# 5:plotting the msa together with the tree#### preparing
# data
t_grouped <- t %>%
    groupOTU(sel_lins, group_name = "SubGsT")
t_grouped@phylo[["tip.label"]] <- sapply(strsplit(t_grouped@phylo[["tip.label"]],
    "\\|"), "[", 1)
p_leg <- ggtree(t_grouped, aes(color = SubGsT), show.legend = F) %<+%
    d + scale_color_manual(values = c("lightgrey", "lightcyan2",
    "bisque3", "darkkhaki")) + geom_nodelab(aes(subset = node ==
    as.numeric(names(sel_lins)[1])), color = "darkblue", nudge_x = -9e-05,
    nudge_y = 30) + geom_nodepoint(aes(subset = node == as.numeric(names(sel_lins)[1])),
    color = "darkblue", shape = 15, size = 4) + geom_nodelab(aes(subset = node ==
    as.numeric(names(sel_lins)[2])), color = "darkred", nudge_x = -9e-05,
    nudge_y = 30) + geom_nodepoint(aes(subset = node == as.numeric(names(sel_lins)[2])),
    color = "darkred", shape = 15, size = 4) + geom_nodelab(aes(subset = node ==
    as.numeric(names(sel_lins)[3])), color = "gold4", nudge_x = -9e-05,
    nudge_y = 30) + geom_nodepoint(aes(subset = node == as.numeric(names(sel_lins)[3])),
    color = "gold4", shape = 15, size = 4) + theme(legend.position = "bottom") +
    guides(fill = guide_legend(nrow = 1)) + new_scale_color()
msap <- msaplot(p_leg, sequences, height = 1, color = c("darkseagreen1",
    "lightblue3", "gray26", "white", "white", "white", "rosybrown2",
    "white")) + scale_fill_manual(values = c("darkseagreen1",
    "lightblue3", "gray26", "white", "white", "white", "rosybrown2",
    "white"), name = "Base", labels = c("A", "C", "G", "N", "N",
    "N", "T", "N")) + theme(legend.title = element_text(face = "bold"))
msapf <- msap + annotate("text", x = seq(0.0013, 0.0013 + 0.000225 *
    5, 0.000225), y = length(t_grouped@phylo[["tip.label"]]) +
    70, label = c("3049", "10116", "22298", "23599", "23604",
    "25613"), angle = 45, size = 5, fontface = "bold") + ylim(c(0,
    length(t_grouped@phylo[["tip.label"]]) + 70))
ggsave("msa_tree_Gamma.jpg", plot = msapf, dpi = 300, height = 7.5,
    width = 10)
