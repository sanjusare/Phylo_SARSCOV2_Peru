# 1:set working directory####
args <- commandArgs(trailingOnly = TRUE)
# 2:load libraries and establish settings####
library(ggplot2)
library(ggtree)
library(treeio)
library(ggpubr)
library(phytools)
library(dplyr)
library(ggthemes)
theme_set(theme_bw(base_size = 15))
Sys.setlocale("LC_TIME", "English")
colors <- c("#cbcbf7", "#f7cbcb", "#f7f4cb", "thistle1", "#91e3db",
    "#abc4d9", "#9c9181", "aquamarine4", "#a9c7af", "#99a5d6",
    "#f4b771", "#ec74dc", "#9ed5ff", "#dedede")
count_colors <- c(Argentina = "dodgerblue", Chile = "forestgreen",
    Peru = "firebrick", Colombia = "orangered", Ecuador = "goldenrod",
    Mexico = "darkmagenta")
# 3:read files#### read maximum clade credibility trees
files <- args[1]
# read trees from tree samples
filetrees <- args[2]
trees_list <- list()
trees_text_list <- list()
meta_list <- list()
for (i in 1:length(filetrees)) {
    name <- as.character(sapply(strsplit(filetrees[i], "_"),
        "[", 1))
    print(name)
    trees_list[[name]] <- read.beast(filetrees[i])
    trees_text_list[[name]] <- readLines(filetrees[i])
    meta_list[[name]] <- read.csv(paste("sample_v2_", name, ".tsv",
        sep = ""), sep = "\t")
}
# 4:Root location probability plot####
height_list <- list()
prob_list <- list()
for (i in 1:length(files)) {
    tr <- read.beast(files[i])
    data <- data.frame(tip = tr@phylo[["tip.label"]], date = as.Date(sapply(strsplit(tr@phylo[["tip.label"]],
        "\\|"), "[", 4)))
    mrt <- Date2decimal(max(data$date))
    locset <- tr@data$location.set[tr@data$node == length(tr@phylo[["tip.label"]]) +
        1][[1]]
    locprob <- tr@data$location.set.prob[tr@data$node == length(tr@phylo[["tip.label"]]) +
        1][[1]]
    height_medi <- as.numeric(tr@data$CAheight_median[tr@data$node ==
        length(tr@phylo[["tip.label"]]) + 1][[1]])
    height_mean <- as.numeric(tr@data$CAheight_mean[tr@data$node ==
        length(tr@phylo[["tip.label"]]) + 1][[1]])
    height_low95 <- as.numeric(tr@data$CAheight_0.95_HPD[tr@data$node ==
        length(tr@phylo[["tip.label"]]) + 1][[1]][1])
    height_hig95 <- as.numeric(tr@data$CAheight_0.95_HPD[tr@data$node ==
        length(tr@phylo[["tip.label"]]) + 1][[1]][2])
    date_medi <- decimal2Date(mrt - height_medi)
    date_mean <- decimal2Date(mrt - height_mean)
    date_low95 <- decimal2Date(mrt - height_low95)
    date_hig95 <- decimal2Date(mrt - height_hig95)
    height_list[[as.character(i)]] <- data.frame(height_median = height_medi,
        height_mean = height_mean, height_low95 = height_low95,
        height_hig95 = height_hig95, date_median = date_medi,
        date_mean = date_mean, date_low95 = date_low95, date_hig95 = date_hig95,
        sample = as.character(sapply(strsplit(files[i], "_"),
            "[", 1)), mrt = mrt)
    prob_list[[as.character(i)]] <- data.frame(country = locset,
        prob = locprob, sample = as.character(sapply(strsplit(files[i],
            "_"), "[", 1)))
}
heights <- do.call(rbind, height_list)
probs <- do.call(rbind, prob_list)
pprobs <- ggplot(probs, aes(x = as.factor(sample), y = prob,
    fill = country)) + geom_bar(stat = "identity", position = "dodge") +
    labs(x = "Sample", y = "Location probability of root") +
    scale_fill_manual(values = c(Argentina = "dodgerblue", Chile = "forestgreen",
        Peru = "firebrick")) + scale_x_discrete(labels = c(1)) +
    theme(legend.title = element_blank(), legend.position = c(0.5,
        0.85)) + ylim(c(0, 1)) + guides(fill = guide_legend(ncol = 3,
    nrow = 1))
pheights <- ggplot(heights, aes(x = as.factor(sample), y = as.POSIXct(date_median))) +
    geom_pointrange(aes(ymin = as.POSIXct(date_low95), ymax = as.POSIXct(date_hig95)),
        color = "firebrick") + geom_point(aes(y = as.POSIXct(date_mean)),
    shape = 21, fill = "dodgerblue", size = 5, alpha = 0.5) +
    labs(y = "Root Date") + scale_x_discrete(labels = c(1)) +
    scale_y_datetime(date_breaks = "1 month", date_labels = "%Y-%b") +
    theme(axis.title.x = element_blank())
p_origin <- ggarrange(pheights, pprobs, ncol = 1, nrow = 2, labels = c("A",
    "B"))
ggsave("Date_Location_of_Lambda_origin.jpg", plot = p_origin,
    dpi = 500, height = 7, width = 4)
# 5:write files of the trees depending of the country of
# the root####
for (i in 1:length(trees_text_list)) {
    tree_names_list <- list()
    for (x in 1:length(trees_text_list[[i]])) {
        if (grepl("tree STATE", trees_text_list[[i]][x])) {
            tree_names_list[[as.character(x)]] <- sapply(strsplit(trees_text_list[[i]][x],
                " "), "[", 2)
        }
    }
    dat_trs <- data.frame(line = names(tree_names_list), state = unlist(tree_names_list))
    dat_trs$number <- 1:nrow(dat_trs)
    names(trees_list[[i]]) <- dat_trs$state
    root_per_list <- c()
    root_arg_list <- c()
    root_chi_list <- c()
    for (z in 1:length(trees_list[[i]])) {
        locroot <- trees_list[[i]][[z]]@data$location[trees_list[[i]][[z]]@data$node ==
            as.character(length(trees_list[[i]][[z]]@phylo[["tip.label"]]) +
                1)][[1]]
        if (locroot == "Peru") {
            root_per_list <- append(root_per_list, dat_trs$line[dat_trs$number ==
                z])
        } else if (locroot == "Argentina") {
            root_arg_list <- append(root_arg_list, dat_trs$line[dat_trs$number ==
                z])
        } else {
            root_chi_list <- append(root_chi_list, dat_trs$line[dat_trs$number ==
                z])
        }
    }
    sel_trees_per <- trees_text_list[[i]][-as.numeric(c(root_arg_list,
        root_chi_list))]
    sel_trees_arg <- trees_text_list[[i]][-as.numeric(c(root_per_list,
        root_chi_list))]
    outfile <- paste(names(trees_text_list)[i], "_v2tree_location_f2_",
        sep = "")
    writeLines(sel_trees_per, paste(outfile, "Per.trees", sep = ""))
    writeLines(sel_trees_arg, paste(outfile, "Arg.trees", sep = ""))
}
# 6:generate tables of transitions####
dats_trees_list <- list()
for (y in 1:length(trees_list)) {
    print(y)
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
        dplyr::mutate(sample = names(trees_list)[y])
    dats_un <- dats %>%
        dplyr::filter(loc1 != loc2)
    dats_gr <- dats %>%
        dplyr::group_by(loc1, loc2) %>%
        dplyr::summarise(freq = n()) %>%
        dplyr::filter(loc1 != loc2)
    dats_com <- list(dats, dats_un, dats_gr)
    dats_trees_list[[names(trees_list)[y]]] <- dats_com
}
# 7:generate tables with the first transition from and to
# an specific country in trees with an specific root####
intr_per_ap_list <- list()
intr_per_ac_list <- list()
intr_per_pa_list <- list()
intr_per_pc_list <- list()
intr_per_cp_list <- list()
intr_per_ca_list <- list()
intr_arg_ap_list <- list()
intr_arg_ac_list <- list()
intr_arg_pa_list <- list()
intr_arg_pc_list <- list()
intr_arg_cp_list <- list()
intr_arg_ca_list <- list()
for (x in 1:length(dats_trees_list)) {
    print(x)
    intr_per_ap_list[[x]] <- dats_trees_list[[x]][[2]] %>%
        dplyr::group_by(tree) %>%
        dplyr::filter(locroot == "Peru" & loc1 == "Argentina" &
            loc2 == "Peru") %>%
        dplyr::filter(height2 == max(height2))
    intr_per_ac_list[[x]] <- dats_trees_list[[x]][[2]] %>%
        dplyr::group_by(tree) %>%
        dplyr::filter(locroot == "Peru" & loc1 == "Argentina" &
            loc2 == "Chile") %>%
        dplyr::filter(height2 == max(height2))
    intr_per_pa_list[[x]] <- dats_trees_list[[x]][[2]] %>%
        dplyr::group_by(tree) %>%
        dplyr::filter(locroot == "Peru" & loc1 == "Peru" & loc2 ==
            "Argentina") %>%
        dplyr::filter(height2 == max(height2))
    intr_per_pc_list[[x]] <- dats_trees_list[[x]][[2]] %>%
        dplyr::group_by(tree) %>%
        dplyr::filter(locroot == "Peru" & loc1 == "Peru" & loc2 ==
            "Chile") %>%
        dplyr::filter(height2 == max(height2))
    intr_per_cp_list[[x]] <- dats_trees_list[[x]][[2]] %>%
        dplyr::group_by(tree) %>%
        dplyr::filter(locroot == "Peru" & loc1 == "Chile" & loc2 ==
            "Peru") %>%
        dplyr::filter(height2 == max(height2))
    intr_per_ca_list[[x]] <- dats_trees_list[[x]][[2]] %>%
        dplyr::group_by(tree) %>%
        dplyr::filter(locroot == "Peru" & loc1 == "Chile" & loc2 ==
            "Argentina") %>%
        dplyr::filter(height2 == max(height2))
    intr_arg_ap_list[[x]] <- dats_trees_list[[x]][[2]] %>%
        dplyr::group_by(tree) %>%
        dplyr::filter(locroot == "Argentina" & loc1 == "Argentina" &
            loc2 == "Peru") %>%
        dplyr::filter(height2 == max(height2))
    intr_arg_ac_list[[x]] <- dats_trees_list[[x]][[2]] %>%
        dplyr::group_by(tree) %>%
        dplyr::filter(locroot == "Argentina" & loc1 == "Argentina" &
            loc2 == "Chile") %>%
        dplyr::filter(height2 == max(height2))
    intr_arg_pa_list[[x]] <- dats_trees_list[[x]][[2]] %>%
        dplyr::group_by(tree) %>%
        dplyr::filter(locroot == "Argentina" & loc1 == "Peru" &
            loc2 == "Argentina") %>%
        dplyr::filter(height2 == max(height2))
    intr_arg_pc_list[[x]] <- dats_trees_list[[x]][[2]] %>%
        dplyr::group_by(tree) %>%
        dplyr::filter(locroot == "Argentina" & loc1 == "Peru" &
            loc2 == "Chile") %>%
        dplyr::filter(height2 == max(height2))
    intr_arg_cp_list[[x]] <- dats_trees_list[[x]][[2]] %>%
        dplyr::group_by(tree) %>%
        dplyr::filter(locroot == "Argentina" & loc1 == "Chile" &
            loc2 == "Peru") %>%
        dplyr::filter(height2 == max(height2))
    intr_arg_ca_list[[x]] <- dats_trees_list[[x]][[2]] %>%
        dplyr::group_by(tree) %>%
        dplyr::filter(locroot == "Argentina" & loc1 == "Chile" &
            loc2 == "Argentina") %>%
        dplyr::filter(height2 == max(height2))
}
intr_per_ap <- do.call(rbind, intr_per_ap_list)
intr_per_ac <- do.call(rbind, intr_per_ac_list)
intr_per_pa <- do.call(rbind, intr_per_pa_list)
intr_per_pc <- do.call(rbind, intr_per_pc_list)
intr_per_cp <- do.call(rbind, intr_per_cp_list)
intr_per_ca <- do.call(rbind, intr_per_ca_list)
intr_arg_ap <- do.call(rbind, intr_arg_ap_list)
intr_arg_ac <- do.call(rbind, intr_arg_ac_list)
intr_arg_pa <- do.call(rbind, intr_arg_pa_list)
intr_arg_pc <- do.call(rbind, intr_arg_pc_list)
intr_arg_cp <- do.call(rbind, intr_arg_cp_list)
intr_arg_ca <- do.call(rbind, intr_arg_ca_list)
intr_1st_all <- do.call(rbind, list(intr_per_ap, intr_per_ac,
    intr_per_pa, intr_per_pc, intr_per_cp, intr_per_ca, intr_arg_ap,
    intr_arg_ac, intr_arg_pa, intr_arg_pc, intr_arg_cp, intr_arg_ca)) %>%
    dplyr::mutate(trans = paste(loc1, "to", loc2))
intr_1st_all$locroot <- factor(intr_1st_all$locroot, levels = c("Peru",
    "Argentina"), ordered = T)
# 8:Plot the distribution of the dates where the first
# transition occurred####
p_trans_pe <- ggplot(subset(intr_1st_all, locroot == "Peru" &
    trans == "Peru to Argentina" | locroot == "Peru" & trans ==
    "Peru to Chile"), aes(x = as.POSIXct(date2), fill = loc2)) +
    geom_histogram(alpha = 0.5, position = position_dodge(),
        bins = 56) + scale_fill_manual(values = count_colors[names(count_colors) %in%
    levels(as.factor(intr_1st_all$loc2))], labels = c("Peru to Argentina",
    "Peru to Chile", "Agentina to Peru")) + labs(y = "Number of trees") +
    theme(legend.position = "none", axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0)) +
    scale_y_continuous(position = "right", limits = c(0, 60)) +
    scale_x_datetime(date_breaks = "1 month", date_labels = "%Y-%b",
        limits = c(as.POSIXct("2020-05-01"), as.POSIXct("2021-02-28"))) +
    geom_vline(data = filter(intr_1st_all, locroot == "Peru"),
        aes(xintercept = as.POSIXct(quantile(intr_1st_all$date2[intr_1st_all$trans ==
            "Peru to Argentina" & intr_1st_all$locroot == "Peru"],
            probs = c(0.05), type = 1)[1])), color = "blue",
        size = 1, alpha = 0.2) + geom_vline(data = filter(intr_1st_all,
    locroot == "Peru"), aes(xintercept = as.POSIXct(quantile(intr_1st_all$date2[intr_1st_all$trans ==
    "Peru to Argentina" & intr_1st_all$locroot == "Peru"], probs = c(0.95),
    type = 1)[1])), color = "blue", size = 1, alpha = 0.2) +
    geom_vline(data = filter(intr_1st_all, locroot == "Peru"),
        aes(xintercept = as.POSIXct(quantile(intr_1st_all$date2[intr_1st_all$trans ==
            "Peru to Chile" & intr_1st_all$locroot == "Peru"],
            probs = c(0.05), type = 1)[1])), color = "darkgreen",
        size = 1, alpha = 0.2) + geom_vline(data = filter(intr_1st_all,
    locroot == "Peru"), aes(xintercept = as.POSIXct(quantile(intr_1st_all$date2[intr_1st_all$trans ==
    "Peru to Chile" & intr_1st_all$locroot == "Peru"], probs = c(0.95),
    type = 1)[1])), color = "darkgreen", size = 1, alpha = 0.2)
p_trans_ar <- ggplot(subset(intr_1st_all, locroot == "Argentina" &
    trans == "Argentina to Peru" | locroot == "Argentina" & trans ==
    "Peru to Chile"), aes(x = as.POSIXct(date2), fill = loc2)) +
    geom_histogram(alpha = 0.5, position = position_dodge(),
        bins = 56) + scale_fill_manual(values = count_colors[names(count_colors) %in%
    levels(as.factor(intr_1st_all$loc2))], labels = c("Peru to Argentina",
    "Peru to Chile", "Agentina to Peru")) + labs(y = "Number of trees") +
    theme(legend.position = "none", axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0)) +
    scale_y_continuous(position = "right", limits = c(0, 60)) +
    scale_x_datetime(date_breaks = "1 month", date_labels = "%Y-%b",
        limits = c(as.POSIXct("2020-05-01"), as.POSIXct("2021-02-28"))) +
    geom_vline(data = filter(intr_1st_all, locroot == "Argentina"),
        aes(xintercept = as.POSIXct(quantile(intr_1st_all$date2[intr_1st_all$trans ==
            "Argentina to Peru" & intr_1st_all$locroot == "Argentina"],
            probs = c(0.05), type = 1)[1])), color = "darkred",
        size = 1, alpha = 0.2) + geom_vline(data = filter(intr_1st_all,
    locroot == "Argentina"), aes(xintercept = as.POSIXct(quantile(intr_1st_all$date2[intr_1st_all$trans ==
    "Argentina to Peru" & intr_1st_all$locroot == "Argentina"],
    probs = c(0.95), type = 1)[1])), color = "darkred", size = 1,
    alpha = 0.2) + geom_vline(data = filter(intr_1st_all, locroot ==
    "Argentina"), aes(xintercept = as.POSIXct(quantile(intr_1st_all$date2[intr_1st_all$trans ==
    "Peru to Chile" & intr_1st_all$locroot == "Argentina"], probs = c(0.05),
    type = 1)[1])), color = "darkgreen", size = 1, alpha = 0.2) +
    geom_vline(data = filter(intr_1st_all, locroot == "Argentina"),
        aes(xintercept = as.POSIXct(quantile(intr_1st_all$date2[intr_1st_all$trans ==
            "Peru to Chile" & intr_1st_all$locroot == "Argentina"],
            probs = c(0.95), type = 1)[1])), color = "darkgreen",
        size = 1, alpha = 0.2)
p_trans <- ggarrange(p_trans_pe, p_trans_ar, ncol = 1, nrow = 2)
ggsave("Transitions_distributions.jpg", plot = p_trans, height = 7.5,
    width = 5)
# 9:generate the map figure####
coord <- read.csv(args[3], sep = "\t", header = FALSE)
names(coord) <- c("Country", "longitude", "latitude")
map_sum <- intr_1st_all %>%
    dplyr::group_by(loc1, loc2, locroot, trans) %>%
    dplyr::summarise(med_date1 = median(date1), med_date2 = median(date2),
        q05date2 = quantile(date2, probs = 0.05, type = 1), q95date2 = quantile(date2,
            probs = 0.95, type = 1)) %>%
    dplyr::left_join(coord, by = c(loc1 = "Country")) %>%
    dplyr::mutate(lat1 = latitude, lon1 = longitude) %>%
    dplyr::select(!c(latitude, longitude)) %>%
    dplyr::left_join(coord, by = c(loc2 = "Country")) %>%
    dplyr::mutate(lat2 = latitude, lon2 = longitude) %>%
    dplyr::select(!c(latitude, longitude))

thismap <- map_data("world")
thismap$fill <- "lightgrey"
thismap$fill[thismap$region == "Peru"] <- "\"firebrick\""
thismap$fill[thismap$region == "Argentina"] <- "\"dodgerblue\""
thismap$fill[thismap$region == "Chile"] <- "\"forestgreen\""

world1 <- ggplot(thismap, aes(long, lat)) + geom_polygon(aes(group = group,
    fill = fill), color = "lightgrey", alpha = 0.5) + scale_fill_manual(values = c("dodgerblue",
    "firebrick", "forestgreen", "grey")) + coord_cartesian(ylim = c(-55,
    0), xlim = c(-80, -55)) + geom_curve(data = subset(map_sum,
    locroot == "Peru" & trans %in% c("Peru to Chile")), aes(x = lat1,
    y = lon1, xend = lat2, yend = lon2), curvature = -0.7, size = 1,
    arrow = arrow(length = unit(0.1, "inches")), color = "darkgreen") +
    geom_curve(data = subset(map_sum, locroot == "Peru" & trans %in%
        c("Peru to Argentina")), aes(x = lat1, y = lon1, xend = lat2,
        yend = lon2), curvature = -0.7, size = 1, arrow = arrow(length = unit(0.1,
        "inches")), color = "blue") + geom_label(data = coord,
    aes(x = latitude, y = longitude, label = Country), nudge_y = -5,
    alpha = 0.5, size = 3, label.padding = unit(0.1, "lines")) +
    theme_map() + theme(legend.position = "none")
world2 <- ggplot(thismap, aes(long, lat)) + geom_polygon(aes(group = group,
    fill = fill), color = "lightgrey", alpha = 0.5) + scale_fill_manual(values = c("dodgerblue",
    "firebrick", "forestgreen", "grey")) + coord_cartesian(ylim = c(-55,
    0), xlim = c(-80, -55)) + geom_curve(data = subset(map_sum,
    locroot == "Argentina" & trans %in% c("Peru to Chile")),
    aes(x = lat1, y = lon1, xend = lat2, yend = lon2), curvature = -0.7,
    size = 1, arrow = arrow(length = unit(0.1, "inches")), color = "darkgreen") +
    geom_curve(data = subset(map_sum, locroot == "Argentina" &
        trans %in% c("Argentina to Peru")), aes(x = lat1, y = lon1,
        xend = lat2, yend = lon2), curvature = 0.7, size = 1,
        arrow = arrow(length = unit(0.1, "inches")), color = "red") +
    geom_label(data = coord, aes(x = latitude, y = longitude,
        label = Country), nudge_y = -5, alpha = 0.5, size = 3,
        label.padding = unit(0.1, "lines")) + theme_map() + theme(legend.position = "none")
p_worlds <- ggarrange(world1, world2, ncol = 1, nrow = 2)
ggsave("Hypotheses_world_maps.jpg", dpi = 300, height = 5, width = 3)
write.table(map_sum, "Summary_of_transitions.tsv", sep = "\t",
    row.names = F, quote = F)
