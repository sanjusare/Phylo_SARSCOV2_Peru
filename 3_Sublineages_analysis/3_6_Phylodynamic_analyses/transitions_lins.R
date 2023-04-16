# 1:set working directory####
args <- commandArgs(trailingOnly = TRUE)
# 2:load libraries and establish settings####
library(ggplot2)
library(ggtree)
library(treeio)
library(ggpubr)
library(phytools)
library(dplyr)
theme_set(theme_bw(base_size = 15))
Sys.setlocale("LC_TIME", "English")
colreg <- c(center = "firebrick", north = "dodgerblue", south = "forestgreen",
    `mid-east` = "darkmagenta", `north-east` = "goldenrod", `south-east` = "hotpink3")
collin <- c(SubL = "forestgreen", SubG1 = "dodgerblue", SubG2 = "firebrick",
    SubG3 = "#ad71ab")
# 3:read trees from the tree samples####
files <- dir(pattern = "*f2.trees")
trees_list <- list()
for (i in 1:length(files)) {
    lin <- as.character(sapply(strsplit(files[i], "_"), "[",
        1))
    subsample <- as.numeric(gsub("v2", "", sapply(strsplit(files[i],
        "_"), "[", 4)))
    name <- paste(lin, subsample, sep = "_")
    print(name)
    trees_list[[name]] <- read.beast(files[i])
}
# 4:generate a table of transitions####
dats_trees_list <- list()
for (y in 1:length(trees_list)) {
    print(names(trees_list)[y])
    dats_list <- list()
    for (z in 1:length(trees_list[[y]])) {
        print(paste(y, z))
        locroot <- trees_list[[y]][[z]]@data$location[trees_list[[y]][[z]]@data$node ==
            as.character(length(trees_list[[y]][[z]]@phylo[["tip.label"]]) +
                1)][[1]]
        heights <- as.data.frame(nodeHeights(trees_list[[y]][[z]]@phylo))
        edges <- as.data.frame(trees_list[[y]][[z]]@phylo[["edge"]])
        names(edges) <- c("node1", "node2")
        names(heights) <- c("height1", "height2")
        heights$height1 <- max(append(heights$height1, heights$height2)) -
            heights$height1
        heights$height2 <- max(append(heights$height1, heights$height2)) -
            heights$height2
        dat <- data.frame(edges, heights)
        maxdate <- Date2decimal(max(as.Date(sapply(strsplit(trees_list[[y]][[z]]@phylo[["tip.label"]],
            "\\|"), "[", 4))))
        dat$he_y1 <- maxdate - heights$height1
        dat$he_y2 <- maxdate - heights$height2
        dat$date1 <- decimal2Date(dat$he_y1)
        dat$date2 <- decimal2Date(dat$he_y2)
        dat$loc1 <- NA
        dat$loc2 <- NA
        trees_list[[y]][[z]]@data$node <- as.numeric(trees_list[[y]][[z]]@data$node)
        for (i in 1:nrow(dat)) {
            dat$loc1[i] <- trees_list[[y]][[z]]@data$location[trees_list[[y]][[z]]@data$node ==
                edges$node1[i]]
            dat$loc2[i] <- trees_list[[y]][[z]]@data$location[trees_list[[y]][[z]]@data$node ==
                edges$node2[i]]
        }
        dat$tree <- z
        dat$locroot <- locroot
        dats_list[[z]] <- dat
    }
    dats <- do.call(rbind, dats_list) %>%
        dplyr::mutate(sample = sapply(strsplit(names(trees_list)[y],
            "_"), "[", 1), subsample = sapply(strsplit(names(trees_list)[y],
            "_"), "[", 2))
    dats_un <- dats %>%
        dplyr::filter(loc1 != loc2)
    dats_gr <- dats %>%
        dplyr::group_by(loc1, loc2) %>%
        dplyr::summarise(freq = n()) %>%
        dplyr::filter(loc1 != loc2) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(perc = freq/sum(freq), sample = sapply(strsplit(names(trees_list)[y],
            "_"), "[", 1), subsample = sapply(strsplit(names(trees_list)[y],
            "_"), "[", 2))
    dats_com <- list(dats, dats_un, dats_gr)
    dats_trees_list[[names(trees_list)[y]]] <- dats_com
}
# 5:calculation of HPD of number of transitions####
for (df in names(dats_trees_list)) {
    print(df)
    dats_trees_list[[df]][[3]]$HPD_L <- 0
    dats_trees_list[[df]][[3]]$HPD_H <- 0
    dats_trees_list[[df]][[3]]$Median <- 0
    dats_trees_list[[df]][[3]]$from_to <- paste(dats_trees_list[[df]][[3]]$loc1,
        "_to_", dats_trees_list[[df]][[3]]$loc2, sep = "")
    test <- dats_trees_list[[df]][[2]] %>%
        dplyr::group_by(tree, loc1, loc2) %>%
        dplyr::summarise(transitions = n()) %>%
        dplyr::ungroup() %>%
        tidyr::complete(tree, loc1, loc2, fill = list(transitions = 0)) %>%
        dplyr::filter(loc1 != loc2) %>%
        dplyr::mutate(from_to = paste(loc1, "_to_", loc2, sep = ""))
    for (i in levels(as.factor(test$from_to))) {
        print(i)
        d <- test %>%
            dplyr::filter(from_to == i)
        dats_trees_list[[df]][[3]]$HPD_L[dats_trees_list[[df]][[3]]$from_to ==
            i] <- round(quantile(d$transitions, 0.025), 0)
        dats_trees_list[[df]][[3]]$HPD_H[dats_trees_list[[df]][[3]]$from_to ==
            i] <- round(quantile(d$transitions, 0.975), 0)
        dats_trees_list[[df]][[3]]$Median[dats_trees_list[[df]][[3]]$from_to ==
            i] <- median(d$transitions)
    }
    all<-dats_trees_list[[df]][[2]]%>%
      dplyr::group_by(tree)%>%
      dplyr::summarise(transitions=n())%>%
      dplyr::ungroup()
    dats_trees_list[[df]][[3]]$HPD_L_tot <- round(quantile(all$transitions, 0.025), 0)
    dats_trees_list[[df]][[3]]$HPD_H_tot <- round(quantile(all$transitions, 0.975), 0)
    dats_trees_list[[df]][[3]]$Median_tot <- median(all$transitions)
}
# 6:heatmap of transitions by region####
medians <- c()
for (i in 1:length(dats_trees_list)) {
    medians <- append(medians, dats_trees_list[[i]][[3]]$Median)
}
p_heat_leg <- ggplot(dats_trees_list[[i]][[3]], aes(x = loc2,
    y = loc1, fill = Median)) + geom_tile() + geom_text(aes(label = paste(HPD_L,
    "-", HPD_H, sep = ""))) + scale_fill_gradient2(low = "darkseagreen",
    high = "cornflowerblue", mid = "white", midpoint = (max(medians) +
        min(medians))/2, limit = c(min(medians), max(medians))) +
    theme(panel.grid = element_blank(), panel.background = element_rect(fill = "lightgrey"),
        axis.title = element_blank(), axis.text.x = element_text(angle = 90,
            hjust = 0, vjust = 0)) + ggtitle(names(dats_trees_list)[i])
heat_leg <- get_legend(p_heat_leg)
heat_ps <- list()
for (i in 1:length(dats_trees_list)) {
    name1<-paste(paste("SubL1+L2", sapply(strsplit(names(dats_trees_list)[i], "_"), "[", 2), sep = "_"), "/",
               dats_trees_list[[i]][[3]]$Median_tot[[1]], "/",
               dats_trees_list[[i]][[3]]$HPD_L_tot[[1]], "-",
               dats_trees_list[[i]][[3]]$HPD_H_tot[[1]], sep = "")
    name2<-paste(names(dats_trees_list)[i], "/",
               dats_trees_list[[i]][[3]]$Median_tot[[1]], "/",
               dats_trees_list[[i]][[3]]$HPD_L_tot[[1]], "-",
               dats_trees_list[[i]][[3]]$HPD_H_tot[[1]], sep = "")
    name<-ifelse(grepl("SubL", names(dats_trees_list)[i]), name1, name2)
    heat_ps[[names(dats_trees_list)[i]]] <- ggplot(dats_trees_list[[i]][[3]],
        aes(x = loc2, y = loc1, fill = Median)) + geom_tile() +
        geom_text(aes(label = paste(HPD_L, "-", HPD_H, sep = ""))) +
        scale_fill_gradient2(low = "darkseagreen", high = "cornflowerblue",
            mid = "white", midpoint = (max(medians) + min(medians))/2,
            limit = c(min(medians), max(medians))) + theme(legend.position = "none",
        panel.grid = element_blank(), panel.background = element_rect(fill = "lightgrey"),
        axis.title = element_blank(), axis.text.x = element_text(angle = 90,
            hjust = 0, vjust = 0)) +
        ggtitle(name)
}
fin_p <- ggarrange(plotlist = c(heat_ps[4], heat_ps[1:3]), ncol = 2,
    nrow = 2, legend.grob = heat_leg, legend = "right") %>%
    annotate_figure(left = text_grob("From", rot = 90, face = "bold",
        size = 20), bottom = text_grob("To", face = "bold", size = 20))
ggsave("Matrices_of_transitions_by_sublineages.jpg", plot = fin_p,
    dpi = 500, height = 7, width = 8)
