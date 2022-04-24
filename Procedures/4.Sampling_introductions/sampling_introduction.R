args<-commandArgs(trailingOnly = TRUE)
Sys.setlocale("LC_TIME", "English")
if(length(args)>6| length(args)<6){
  print("usage=Rscript sampling_introductions.R [metadata_file] [seqs_list_file] [country] [variant] [days_before] [days_after]")
  break
}
if (!require(dplyr, quietly = T)) install.packages("dplyr")
library(dplyr)
if (!require(tidyr, quietly = T)) install.packages("tidyr")
library(tidyr)
tips<-readLines(args[2])
metadata<-read.csv(file=args[1],sep="\t",header=TRUE)%>%
  dplyr::filter(name_analysis%in%tips)
metadata$date<-as.Date(metadata$date)
take_lin<-function(lin,begin,end,location){
  lin_peru<-metadata%>%
    dplyr::filter(grepl(lin,variant)&Country==location)%>%
    dplyr::select(name_analysis,pango_lineage,variant,Country,date)%>%
    dplyr::mutate(date=as.Date(date))%>%
    dplyr::arrange(date)
  lin_tsv<-metadata%>%
    dplyr::filter(grepl(lin,variant))%>%
    dplyr::filter(date>=(lin_peru$date[1]-begin))%>%
    dplyr::filter(date<=(lin_peru$date[1]+end))%>%
    dplyr::select(name_analysis,pango_lineage,variant,Country,date)
  lin_tsv_u<-lin_tsv%>%
    distinct(name_analysis,.keep_all=TRUE)
  #lin_tsv_u tiene todos los genomas de la variante de interes desde "begin" dias antes del primer genoma peruano de la variante de interes y "end" genomas despues del primer genoma peruano de la variante de interes
  met_filt_1_unido_meta<-metadata%>%
    filter(name_analysis%in%lin_tsv_u$name_analysis)
  name_list<-paste(gsub("\\.","",lin),"_introduction.list",sep="")
  name_metadata<-paste(gsub("\\.","",lin),"_introduction.tsv",sep="")
  writeLines(met_filt_1_unido_meta$name_analysis,name_list)
  write.table(met_filt_1_unido_meta,name_metadata,quote=FALSE,sep="\t",col.names=TRUE,row.names=FALSE)
  return(met_filt_1_unido_meta)
}
Gamma<-take_lin(args[4],args[5],args[6],args[3])
