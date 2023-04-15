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
colreg <- c(center = "firebrick", north = "dodgerblue", south = "forestgreen",
    `mid-east` = "darkmagenta", `north-east` = "goldenrod", `south-east` = "hotpink3")
collin <- c(SubG1 = "dodgerblue", SubG2 = "firebrick", SubG3 = "#ad71ab")
# 3:read files#### read the maximum likelihoood tree
t <- read.iqtree(args[1])
# read the metadata of genomes from Peru and filter those
# that are present in the tree
d <- read.csv(args[2], sep = "\t") %>%
    dplyr::filter(name_analysis %in% t@phylo[["tip.label"]]) %>%
    mutate(date = as.Date(date))
# read the file of the grouping information of cities in
# regions
reg <- read.csv(args[3], sep = "\t")
# 4:Sublineages identification#### metadata arrangements
d$City[d$City == "Amazona"] <- "Amazonas"
d$City[d$City == "HuancavelIca"] <- "Huancavelica"
d$City[d$City == "Huánuco"] <- "Huanuco"
d$City[d$City == "Lalibertad"] <- "LaLibertad"
d$City[d$City == "MadreDeDios"] <- "MadredeDios"
d <- dplyr::left_join(d, reg)
# selection of sublineages
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
# save the tips that belong to each sublineage
for (i in 1:length(sel_lins)) {
    name_lin <- paste("sublineage_SubG", i, ".list", sep = "")
    writeLines(sel_lins[[i]], name_lin)
}
# 5:plotting the sublineages in the tree#### preparing data
t_grouped <- t %>%
    groupOTU(sel_lins, group_name = "SubGsT")
p_leg <- ggtree(t_grouped, aes(color = SubGsT), show.legend = F) %<+%
    d + scale_color_manual(values = c("lightgrey", "lightcyan2",
    "bisque3", "darkkhaki")) + new_scale_color() + geom_tippoint(aes(subset = label %in%
    unlist(sel_lins), color = Region)) + scale_color_manual(values = colreg) +
    geom_nodelab(aes(subset = node == as.numeric(names(sel_lins)[1])),
        color = "darkblue", nudge_x = -9e-05, nudge_y = 20) +
    geom_nodepoint(aes(subset = node == as.numeric(names(sel_lins)[1])),
        color = "darkblue", shape = 15, size = 4) + geom_nodelab(aes(subset = node ==
    as.numeric(names(sel_lins)[2])), color = "darkred", nudge_x = -9e-05,
    nudge_y = 20) + geom_nodepoint(aes(subset = node == as.numeric(names(sel_lins)[2])),
    color = "darkred", shape = 15, size = 4) + geom_nodelab(aes(subset = node ==
    as.numeric(names(sel_lins)[3])), color = "gold4", nudge_x = -9e-05,
    nudge_y = 20) + geom_nodepoint(aes(subset = node == as.numeric(names(sel_lins)[3])),
    color = "gold4", shape = 15, size = 4) + theme(legend.position = "bottom") +
    guides(color = guide_legend(nrow = 1))
leg <- get_legend(p_leg)
# specify that the data is from Gamma
t_grouped_gam <- t_grouped
d_gam <- d
sel_lins_gam <- sel_lins
# plotting the tree
p_t_fin_gam <- p_leg + theme(legend.position = "none")
ggsave("Maximum_likelihood_Gamma.jpg", plot = p_t_fin_gam, dpi = 500,
    height = 5, width = 5)
