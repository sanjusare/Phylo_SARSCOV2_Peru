# 1:set working directory####
args <- commandArgs(trailingOnly = TRUE)
# 2:load libraries and establish settings####
library(dplyr)
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
metadata <- read.csv(args[2], sep = "\t", header = TRUE)
# 4:create the list of Lambda genomes from the selected
# countries to extract from GISAID####
datos_sel <- metadata %>%
    filter(grepl("Lambda", variant)) %>%
    filter(Country %in% c("Argentina", "Chile", "Peru", "Ecuador",
        "Colombia", "Mexico")) %>%
    mutate(date = as.Date(date)) %>%
    dplyr::filter(date >= as.Date("2021-01-15")) %>%
    dplyr::filter(date <= as.Date("2021-12-31")) %>%
    dplyr::select(name_analysis, pango_lineage, variant, Country,
        date) %>%
    dplyr::mutate(week = as.Date(cut.Date(date, breaks = "1 week")))
met_filt <- metadata %>%
    filter(name_analysis %in% datos_sel$name_analysis)
met_filt$age <- gsub("#N/A", "Unknown", met_filt$age)
tabla_rename <- met_filt %>%
    dplyr::select(strain, name_analysis)
write.table(met_filt, "sample_initial.tsv", quote = FALSE, sep = "\t",
    col.names = TRUE, row.names = FALSE)
writeLines(met_filt$strain, "sample_initial.list")
write.table(tabla_rename, "sample_rename_initial.tsv", sep = "\t",
    quote = FALSE, col.names = FALSE, row.names = FALSE)

