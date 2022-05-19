args<-commandArgs(trailingOnly = TRUE)
Sys.setlocale("LC_TIME", "English")
if(length(args)>1| length(args)<1){
  print("usage=nexclean.R [nexclade_file]")
  break
}
if (!require(dplyr, quietly = T)) install.packages("dplyr")
library(dplyr)
pat<-paste(args[1],"*",sep="")
tabs<-dir(pattern=pat)
tabs_list<-list()
for(i in 1:length(tabs)){
  tabs_list[[i]]<-read.csv(tabs[i],sep=";")
}
tab<-do.call(rbind,tabs_list)
tab_clean<-tab%>%
  dplyr::filter(qc.overallStatus=="good")
tab_clean$VOC<-tab_clean$clade
for(i in 1:length(tab_clean$clade)){
  if(grepl("Omicron",tab_clean$clade[i])){
    tab_clean$clade[i]<-"Omicron"
  }else if(grepl("Delta",tab_clean$clade[i])){
    tab_clean$clade[i]<-"Delta"
  }else if(grepl("Alpha",tab_clean$clade[i])){
    tab_clean$clade[i]<-"Alpha"
  }else if(grepl("Mu",tab_clean$clade[i])){
    tab_clean$clade[i]<-"Mu"
  }else if(grepl("Lambda",tab_clean$clade[i])){
    tab_clean$clade[i]<-"Lambda"
  }else if(grepl("Gamma",tab_clean$clade[i])){
    tab_clean$clade[i]<-"Gamma"
  }else{
    tab_clean$clade[i]<-"Other"
  }
}
tab_final<-tab_clean%>%
  dplyr::filter(clade!="Other")
writeLines(tab_final$seqName,"nexclean.list")