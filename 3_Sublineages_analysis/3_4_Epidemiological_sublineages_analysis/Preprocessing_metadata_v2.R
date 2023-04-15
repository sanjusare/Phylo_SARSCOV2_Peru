args <- commandArgs(trailingOnly = TRUE)
library(dplyr)
library(tidyr)
tabla <- read.csv(args[1], sep = "\t")
tabla$strain <- paste(tabla$Virus.name, tabla$Collection.date,
    tabla$Submission.date, sep = "|")
tabla$Location <- gsub(" ", "", tabla$Location)
tabla$Continent <- sapply(strsplit(tabla$Location, "/"), "[",
    1)
tabla$Country <- sapply(strsplit(tabla$Location, "/"), "[", 2)
tabla$City <- sapply(strsplit(tabla$Location, "/"), "[", 3)
tabla$Collection_date <- as.Date(tabla$Collection.date)
tabla$name_analysis <- paste(tabla$Accession.ID, tabla$Continent,
    tabla$Country, tabla$Collection_date, sep = "|")
tabla_u <- tabla %>%
    dplyr::distinct(strain, .keep_all = TRUE) %>%
    tidyr::drop_na(Collection_date) %>%
    dplyr::filter(Is.complete. == "True") %>%
    dplyr::filter(Host == "Human")
final_table <- tabla_u %>%
    dplyr::select(name_analysis, strain, EPI_ISL = Accession.ID,
        date = Collection_date, host = Host, age = Patient.age,
        gender = Gender, pango_lineage = Pango.lineage, pango_ver = Pangolin.version,
        variant = Variant, Continent, Country, City)
write.table(final_table, "ordered_metadata.tsv", sep = "\t",
    quote = FALSE, row.names = FALSE)
