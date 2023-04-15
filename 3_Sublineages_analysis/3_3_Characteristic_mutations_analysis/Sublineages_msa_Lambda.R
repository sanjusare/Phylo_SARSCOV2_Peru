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
collin <- c(SubL1 = "#58578a", SubL2 = "#8ad3db")
# 3:read files####
args <- c("Lambda_good_seqs.fasta.treefile", "nexclean_Peru_Var.tsv",
    "Lambda_good_seqs_spec_sites.fasta")
# read the maximum likelihoood tree
t <- read.iqtree(args[1])
# read the metadata of genomes from Peru and filter those
# that are present in the tree
d <- read.csv(args[2], sep = "\t") %>%
    dplyr::filter(name_analysis %in% t@phylo[["tip.label"]]) %>%
    mutate(date = as.Date(date))
# read the fasta of the selected positions
sequences <- read.fasta(args[3])
for (i in 1:length(sequences)) {
    sequences[[i]][2] <- sequences[["EPI_ISL_5934945"]][1]
    sequences[[i]][3] <- sequences[["EPI_ISL_5934945"]][1]
}
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
    groupOTU(sel_lins, group_name = "SubLsT")
t_grouped@phylo[["tip.label"]] <- sapply(strsplit(t_grouped@phylo[["tip.label"]],
    "\\|"), "[", 1)
p_leg <- ggtree(t_grouped, aes(color = SubLsT), show.legend = F) %<+%
    d + scale_color_manual(values = c("lightgrey", "#58578a",
    "#8ad3db")) + geom_nodelab(aes(subset = node == as.numeric(names(sel_lins)[1])),
    color = "#58578a", nudge_x = -8e-05, nudge_y = 40) + geom_nodepoint(aes(subset = node ==
    as.numeric(names(sel_lins)[1])), color = "#58578a", shape = 15,
    size = 4) + geom_nodelab(aes(subset = node == as.numeric(names(sel_lins)[2])),
    color = "#8ad3db", nudge_x = -8e-05, nudge_y = 40) + geom_nodepoint(aes(subset = node ==
    as.numeric(names(sel_lins)[2])), color = "#8ad3db", shape = 15,
    size = 4) + theme(legend.position = "bottom") + guides(color = guide_legend(nrow = 1)) +
    new_scale_color()
msap <- msaplot(p_leg, sequences, height = 1, color = c("lightblue3",
    "white", "rosybrown2", "white"), window = c(1, 2), width = 0.25) +
    scale_fill_manual(values = c("lightblue3", "white", "rosybrown2",
        "white"), name = "Base", labels = c("C", "", "T", "")) +
    theme(legend.title = element_text(face = "bold"))
msapf <- msap + annotate("text", x = 0.0014, y = length(t_grouped@phylo[["tip.label"]]) +
    150, label = c("28849"), angle = 45, size = 5, fontface = "bold") +
    ylim(c(0, length(t_grouped@phylo[["tip.label"]]) + 150)) +
    xlim(c(0, 0.00175))
ggsave("msa_tree_Lambda.jpg", plot = msapf, dpi = 500, height = 7.5,
    width = 10)
