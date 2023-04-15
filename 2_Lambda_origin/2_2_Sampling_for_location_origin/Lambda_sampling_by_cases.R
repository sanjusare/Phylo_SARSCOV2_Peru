# 1:set working directory####
args <- commandArgs(trailingOnly = TRUE)
# 2:load libraries and establish settings####
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(ggpmisc)
theme_set(theme_bw(base_size = 15))
Sys.setlocale("LC_TIME", "English")
cols <- c(Omicron = "#db8a8a", Gamma = "#dbd48a", Lambda = "#8fdb8a",
    Delta = "#8ad3db", Mu = "#a18adb", Alpha = "#db8ac0", Other = "darkgrey")
colors <- c("#cbcbf7", "#f7cbcb", "#f7f4cb", "thistle1", "#91e3db",
    "#abc4d9", "#9c9181", "aquamarine4", "#a9c7af", "#99a5d6",
    "#f4b771", "#ec74dc", "#9ed5ff", "#dedede")
count_colors <- c(Argentina = "dodgerblue", Chile = "forestgreen",
    Peru = "firebrick", Colombia = "orangered", Ecuador = "goldenrod",
    Mexico = "darkmagenta")
# 3:read files####
data <- read.csv(file = args[1], sep = "\t", header = TRUE)
data$date <- as.Date(data$date)
final_data <- read.csv(file = args[2], sep = "\t", header = TRUE)
# 4:Sampling function####
sampling_2 <- function(seed = 9876345, number = "19", from = "2020-12-01",
    to = "2021-09-30", genpcas = 1000) {
    samp_dat <- final_data %>%
        dplyr::mutate(collection_date = as.Date(collection_date)) %>%
        dplyr::filter(collection_date >= as.Date(from) & collection_date <=
            as.Date(to) & Country %in% c("Argentina", "Chile",
            "Peru")) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(samp_size = round(Lambda_cases/(sum(Lambda_cases,
            na.rm = T)/genpcas), 0))
    set.seed(seed)
    datos_tsv <- data %>%
        dplyr::filter(date >= as.Date(from)) %>%
        dplyr::filter(date <= as.Date(to)) %>%
        dplyr::select(name_analysis, pango_lineage, variant,
            Country, date) %>%
        filter(Country %in% c("Argentina", "Chile", "Peru")) %>%
        group_by(Country) %>%
        dplyr::mutate(week = as.Date(cut.Date(date, breaks = "1 week")))
    datos_tsv_sum <- datos_tsv %>%
        dplyr::group_by(week, Country) %>%
        dplyr::summarise(new_count = n())
    datos_list <- list()
    for (i in 1:nrow(samp_dat)) {
        print(i)
        name <- paste(samp_dat$Country[i], samp_dat$collection_date[i],
            sep = "_")
        if (samp_dat$samp_size[i] > 0) {
            datos <- datos_tsv %>%
                dplyr::filter(Country == samp_dat$Country[i] &
                  week == samp_dat$collection_date[i])
            if (nrow(datos) >= samp_dat$samp_size[i]) {
                datos_list[[name]] <- datos %>%
                  dplyr::sample_n(samp_dat$samp_size[i])
            } else {
                datos_list[[name]] <- datos
            }
        }
    }
    datos_fin <- do.call(rbind, datos_list)
    met_filt <- data %>%
        filter(name_analysis %in% datos_fin$name_analysis)
    met_filt$age <- gsub("#N/A", "Unknown", met_filt$age)
    write.table(met_filt, paste("sample_v2_", number, ".tsv",
        sep = ""), quote = FALSE, sep = "\t", col.names = TRUE,
        row.names = FALSE)
    writeLines(met_filt$name_analysis, paste("sample_v2_", number,
        ".list", sep = ""))
    gen_samp <- datos_fin %>%
        group_by(week, Country) %>%
        summarise(samp_gen = n()) %>%
        right_join(samp_dat, by = c(week = "collection_date",
            Country = "Country")) %>%
        left_join(datos_tsv_sum)
    gen_samp$samp_gen[is.na(gen_samp$samp_gen)] <- 0
    factor <- max(samp_dat$count, na.rm = T)
    p <- ggplot(gen_samp, aes(x = as.POSIXct(week))) + geom_bar(aes(y = new_count),
        stat = "identity", fill = "lightgrey") + geom_bar(aes(y = samp_gen,
        fill = Country), stat = "identity") + geom_bar(data = subset(gen_samp,
        samp_size > 0), aes(y = samp_size), stat = "identity",
        color = "red", fill = "transparent", linewidth = 0.2) +
        geom_line(aes(y = Lambda_cases/factor)) + facet_wrap(~Country,
        nrow = 3, ncol = 1, scales = "free_y") + scale_y_continuous(sec.axis = sec_axis(trans = ~. *
        factor, name = "Lambda cases")) + scale_x_datetime(date_breaks = "1 month",
        date_labels = "%Y-%b", limits = c(as.POSIXct(from), as.POSIXct(to))) +
        theme(legend.position = "none", axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 0,
                vjust = 0.5)) + scale_fill_manual(values = count_colors)
    ggsave(paste("sample_v2_", number, ".jpg", sep = ""), plot = p,
        dpi = 300, height = 9, width = 5)
    return(gen_samp)
}
d <- sampling_2(args[3], args[4], args[5], args[6], args[7])
# 5:Plot of correlations#### before sampling
cor <- ggplot(d, aes(x = Lambda_cases, y = new_count)) + geom_point(aes(color = Country)) +
    scale_color_manual(values = count_colors[names(count_colors) %in%
        levels(as.factor(d$Country))]) + geom_smooth(method = "lm") +
    stat_poly_eq(aes(label = after_stat(rr.label))) + labs(y = "Lambda Genomes",
    x = "Lambda cases") + theme(legend.position = "bottom", legend.title = element_blank()) +
    guides(color = guide_legend(nrow = 1))
leg <- get_legend(cor)
cor_notleg <- cor + theme(legend.position = "none") + ggtitle("Before sampling")
# after sampling
cor_sam <- ggplot(d, aes(x = Lambda_cases, y = samp_gen)) + geom_point(aes(color = Country)) +
    scale_color_manual(values = count_colors[names(count_colors) %in%
        levels(as.factor(d$Country))]) + geom_smooth(method = "lm") +
    stat_poly_eq(aes(label = after_stat(rr.label))) + labs(y = "Lambda Genomes",
    x = "Lambda cases") + theme(legend.position = "none") + ggtitle("After sampling")
# final cor plot
p_cors <- ggarrange(cornotleg, corsam, legend = "bottom", legend.grob = leg,
    common.legend = T)
ggsave("Lambda_correlation_plot.jpg", plot = p_cors, dpi = 500,
    height = 7.5, width = 5)
