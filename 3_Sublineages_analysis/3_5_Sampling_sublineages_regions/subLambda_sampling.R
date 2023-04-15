# 1:set working directory####
args <- commandArgs(trailingOnly = TRUE)
# 2:load libraries and establish settings####
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
theme_set(theme_bw(base_size = 12))
Sys.setlocale("LC_TIME", "English")
colreg <- c(center = "firebrick", north = "dodgerblue", south = "forestgreen",
    `mid-east` = "darkmagenta", `north-east` = "goldenrod", `south-east` = "hotpink3")
# 3:read files#### read the epidemiological data (from
# section 3_4, step 5)
f_df_cas_all <- read.csv(file = args[1], sep = "\t", header = TRUE) %>%
    dplyr::mutate(collection_date = as.Date(collection_date))
# read the genome data with sublineages information (from
# section 3_4, step 5)
df_p_reg <- read.csv(file = args[2], sep = "\t", header = TRUE)
# 4:joining the data of the two sublineages of Lambda####
f_df_cas_all_lam <- f_df_cas_all %>%
    dplyr::ungroup() %>%
    dplyr::group_by(week, Region) %>%
    dplyr::mutate(lam_cases = sum(smooth_lin_cases), lam_count = sum(count)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(SubLs == "SubL1") %>%
    dplyr::mutate(samp_gen = round(lam_cases/(sum(lam_cases,
        na.rm = T)/200), 0), collection_date = as.Date(collection_date))
f_df_cas_all_lam$samp_gen[is.na(f_df_cas_all_lam$samp_gen)] <- 0
# 5:sampling the data according to cases by region and
# week####
datos_list <- list()
for (z in 1:nrow(f_df_cas_all_lam)) {
    print(z)
    set.seed(234364151)
    if (f_df_cas_all_lam$samp_gen[z] > 0) {
        name <- paste(f_df_cas_all_lam$Region[z], f_df_cas_all_lam$collection_date[z],
            sep = "_")
        te <- df_p_reg %>%
            dplyr::mutate(collection_date = as.Date(collection_date)) %>%
            dplyr::filter(Region == f_df_cas_all_lam$Region[z] &
                collection_date == f_df_cas_all_lam$collection_date[z] &
                SubLs %in% c("SubL1", "SubL2"))
        if (nrow(te) >= f_df_cas_all_lam$samp_gen[z]) {
            datos_list[[name]] <- te %>%
                dplyr::sample_n(f_df_cas_all_lam$samp_gen[z])
        } else {
            datos_list[[name]] <- te
        }
    }
}
datos_lam <- do.call(rbind, datos_list)
gen_group_cases <- datos_lam %>%
    dplyr::group_by(Region, collection_date) %>%
    dplyr::summarise(samp_count = n()) %>%
    dplyr::right_join(f_df_cas_all_lam)
gen_group_cases$samp_count[is.na(gen_group_cases$samp_count)] <- 0
# 6:plot of correlation between estimated cases and number
# of genomes by week by region before and after
# sampling####
p_cor_lam <- ggscatter(subset(gen_group_cases, samp_count !=
    0), x = "lam_count", y = "lam_cases", color = "black", fill = "Region",
    shape = 21, size = 3, add = "reg.line", add.params = list(color = "blue",
        fill = "lightgray"), conf.int = TRUE, cor.coef = TRUE,
    cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")) +
    labs(y = "Number of Cases", x = "Number of genomes") + scale_fill_manual(values = colreg) +
    theme(legend.position = "bottom") + guides(fill = guide_legend(nrow = 1))
cor_leg <- get_legend(p_cor_lam)
p_cor_l <- p_cor_lam + theme(legend.position = "none")
p_cor_sam_lam <- ggscatter(subset(gen_group_cases, samp_count !=
    0), x = "samp_count", y = "lam_cases", color = "black", fill = "Region",
    shape = 21, size = 3, add = "reg.line", add.params = list(color = "blue",
        fill = "lightgray"), conf.int = TRUE, cor.coef = TRUE,
    cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")) +
    labs(y = "Number of Cases", x = "Number of genomes") + scale_fill_manual(values = colreg) +
    theme(legend.position = "none") + guides(fill = guide_legend(nrow = 1))
p_cor_fin <- ggarrange(p_cor_lam, p_cor_sam_lam, labels = c("A",
    "B"), nrow = 2, ncol = 1, common.legend = T, legend = "bottom",
    legend.grob = cor_leg)
ggsave("Lambda_sampling_correlation.jpg", plot = p_cor_fin, dpi = 500,
    height = 10, width = 7.5)
# 7:Plot to evaluate cases, available genomes and sampled
# genomes####
factor <- max(gen_group_cases$lam_count, na.rm = T)
p_sampled <- ggplot(gen_group_cases, aes(x = as.POSIXct(collection_date))) +
    geom_bar(aes(y = lam_count), stat = "identity", fill = "lightgrey") +
    geom_bar(aes(y = samp_count, fill = Region), stat = "identity") +
    geom_bar(data = subset(gen_group_cases, samp_gen > 0), aes(y = samp_gen),
        stat = "identity", color = "red", fill = "transparent",
        linewidth = 0.2) + geom_line(aes(y = lam_cases/factor)) +
    facet_wrap(~Region, nrow = 3, ncol = 2, scales = "free_y") +
    scale_y_continuous(sec.axis = sec_axis(trans = ~. * factor,
        name = "SubL1+L2 cases")) + scale_x_datetime(date_breaks = "1 month",
    date_labels = "%Y-%b", limits = c(as.POSIXct("2020-10-01"),
        as.POSIXct("2022-01-31"))) + theme(legend.position = "none",
    axis.title.x = element_blank(), axis.text.x = element_text(angle = 90,
        hjust = 0, vjust = 0.5)) + scale_fill_manual(values = colreg) +
    labs(y = "Number of genomes")
ggsave("Lambda_sampling_distribution.jpg", plot = p_sampled,
    dpi = 300, height = 7.5, width = 7.5)
# 8:Saving the sampled genomes and tables with relevant
# information for analysis####
d_reg <- datos_lam %>%
    dplyr::select(name_analysis, Region)
writeLines(datos_lam$name_analysis, paste("Lambda", "_BI_fin_v2.list",
    sep = ""))
write.table(datos_lam, paste("Lambda", "_meta_fin_v2.tsv", sep = ""),
    sep = "\t", quote = F, row.names = F)
write.table(d_reg, paste("Lambda", "_region_fin_v2.tsv", sep = ""),
    sep = "\t", quote = F, row.names = F)
