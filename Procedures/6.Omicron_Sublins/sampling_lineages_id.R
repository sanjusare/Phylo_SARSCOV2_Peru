args<-commandArgs(trailingOnly = TRUE)
Sys.setlocale("LC_TIME", "English")
if(length(args)>4| length(args)<4){
  print("usage=Rscript sampling_lineages_id.R [metadata_file] [seqs_list_file] [country] [VOC]")
  break
}
if (!require(dplyr, quietly = T)) install.packages("dplyr")
library(dplyr)
if (!require(tidyr, quietly = T)) install.packages("tidyr")
library(tidyr)
tabla<-read.csv(args[1],sep="\t")
nexlist<-readLines(args[2])
df_p<-tabla%>%
  dplyr::filter(grepl(args[4],variant))%>%
  drop_na(date)%>%
  filter(name_analysis%in%nexlist)
df_p_per<-df_p%>%
  dplyr::filter(Country==args[3])%>%
  dplyr::mutate(date=as.Date(date))%>%
  drop_na(date)
df_p_world<-df_p%>%
  dplyr::filter(Country!=args[3])%>%
  dplyr::mutate(date=as.Date(date))%>%
  drop_na(date)
df_p_world_lin<-df_p_world%>%
  dplyr::group_by(Country,pango_lineage)%>%
  dplyr::summarise(count=n())
df_p_world_lin_list<-list()
for(i in 1:nrow(df_p_world_lin)){
  row<-df_p_world_lin[i,]
  df_p_world_lin_list[[paste(row$Country,row$pango_lineage,sep="_")]]<-df_p_world%>%
    dplyr::filter(Country==row$Country,pango_lineage==row$pango_lineage)%>%
    dplyr::sample_n(1)
}
df_p_world_lin_fin<-do.call(rbind,df_p_world_lin_list)
df_for_iqtree<-rbind(df_p_per,df_p_world_lin_fin)
name<-paste(args[4],"_",args[3],"_","sample_lin_id",sep="")
write.table(df_for_iqtree,paste(name,".tsv",sep=""),sep="\t",col.names=T,row.names=F,quote=F)
writeLines(df_for_iqtree$EPI_ISL,paste(name,".list",sep=""))
