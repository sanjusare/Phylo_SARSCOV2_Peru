args <- commandArgs(trailingOnly = TRUE)
Sys.setlocale("LC_TIME", "English")
if (length(args) > 2 | length(args) < 2) {
    print("usage=all_metadata_fRT.R [metadata] [nextclade_list]")
    break
}
# preparations####
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
Sys.setlocale("LC_TIME", "English")
theme_set(theme_bw(base_size = 12))
cols <- c(Omicron = "#db8a8a", Gamma = "#dbd48a", Lambda = "#8fdb8a",
    Delta = "#8ad3db", Mu = "#a18adb", Alpha = "#db8ac0", "darkgrey")
metadata <- read.csv(file = args[1], sep = "\t", header = TRUE)
nexlist <- readLines(args[2])
f_d_c <- metadata %>%
    dplyr::filter(Country == "Peru") %>%
    mutate(date = as.Date(date)) %>%
    filter(date > "2020-10-01") %>%
    filter(date < "2022-01-01") %>%
    filter(name_analysis %in% nexlist)
write.table(f_d_c, "nexclean_Peru_Var.tsv", sep = "\t", quote = FALSE,
    row.names = FALSE)
