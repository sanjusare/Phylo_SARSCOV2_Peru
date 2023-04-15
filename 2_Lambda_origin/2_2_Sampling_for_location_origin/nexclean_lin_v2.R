args <- commandArgs(trailingOnly = TRUE)
Sys.setlocale("LC_TIME", "English")
if (length(args) > 2 | length(args) < 2) {
    print("usage=nexclean.R [nexclade_file] [metadata_file]")
    break
}
library(dplyr)
library(ggplot2)
name <- sapply(strsplit(args[1], "\\."), "[", 1)
tab <- read.csv(args[1], sep = "\t")
tab_clean <- tab %>%
    dplyr::filter(qc.overallStatus == "good")
writeLines(tab_clean$seqName, paste("nexclean_", name, ".list",
    sep = ""))
data <- read.table(args[2], sep = "\t", quote = "", header = TRUE)
data$name_analysis <- gsub("'", "_", data$name_analysis)
data_clean <- data %>%
    dplyr::filter(name_analysis %in% tab_clean$seqName)
write.table(data_clean, paste("final_metadata_", name, ".tsv",
    sep = ""), quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
data_clean_sel <- data_clean %>%
    select(name_analysis, Country)
data_clean_sel$Country[!(data_clean_sel$Country %in% c("Argentina",
    "Chile", "Peru"))] <- "Other"
write.table(data_clean_sel, paste("location_", name, ".tsv",
    sep = ""), quote = FALSE, sep = "\t", col.names = FALSE,
    row.names = FALSE)
