args<-commandArgs(trailingOnly = TRUE)
library(dplyr)
d<-read.csv(args[1],sep="\t")
d_L<-d%>%
  dplyr::filter(grepl("Lambda",variant))
d_G<-d%>%
  dplyr::filter(grepl("Gamma",variant))
writeLines(d_L$name_analysis,"Lambda_good.list")
writeLines(d_G$name_analysis,"Gamma_good.list")
