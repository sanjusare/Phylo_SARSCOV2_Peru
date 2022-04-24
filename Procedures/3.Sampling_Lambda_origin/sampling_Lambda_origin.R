args<-commandArgs(trailingOnly = TRUE)
Sys.setlocale("LC_TIME", "English")
if(length(args)>4| length(args)<4){
  print("usage=Rscript sampling_Lambda_origin.R [metadata_file] [seqs_list_file] [date_lower_limit] [date_upper_limit]")
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
datos_tsv<-metadata%>%
  dplyr::filter(date>=as.Date(args[3]))%>%
  dplyr::filter(date<=as.Date(args[4]))%>%
  dplyr::select(name_analysis,pango_lineage,variant,Country,date)
#parte1####
datos_tsv_lambda<-datos_tsv%>%
  filter(grepl("Lambda",variant))
#datos_tsv_lambda tiene todos los lambda dentro del rango de fechas seleccionadas
#parte2####
parte_2<-datos_tsv%>%
  filter(!(grepl("Lambda",variant)),Country=="Peru")
#parte_2 contiene todos los genomas peruanos que no son Lambda
#parte3####
pob_1_wor_sum<-datos_tsv%>%
  dplyr::filter(Country != "Peru" & !(grepl("Lambda",variant)))%>%
  group_by(Country,pango_lineage,date)%>%
  drop_na(date)%>%
  summarise(count=n())
pob_1_wor_list<-list()
for (i in 1:nrow(pob_1_wor_sum)){
  print(i)
  row<-pob_1_wor_sum[i,] #crea obj de 1 fila
  name<-paste(row$Country,row$pango_lineage,row$date,sep="_")
  pob_1_wor_list[[name]]<-datos_tsv%>%
    dplyr::filter(Country==row$Country,pango_lineage==row$pango_lineage,date==row$date)%>%
    dplyr::sample_n(size=1)
}
pob_1_wor<-do.call(rbind,pob_1_wor_list)
#pob_1_wor tiene 1 de cada linaje de cada fecha disponible de cada pais diferente a peru diferente a lambda
#parte final####
data_ML1<-rbind(datos_tsv_lambda,parte_2)
data_ML<-rbind(data_ML1,pob_1_wor)
data_ML_u<-data_ML%>%
  distinct(name_analysis,.keep_all=TRUE)
met_filt<-metadata%>%
  filter(name_analysis%in%data_ML_u$name_analysis)
met_filt$age<-gsub("#N/A","Unknown",met_filt$age)
write.table(met_filt,"Lambda_origin.tsv",quote=FALSE,sep="\t",col.names=TRUE,row.names=FALSE)
writeLines(met_filt$name_analysis,"Lambda_origin.list")