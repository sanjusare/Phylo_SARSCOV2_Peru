args<-commandArgs(trailingOnly = TRUE)
Sys.setlocale("LC_TIME")
if(length(args) != 8){
  print("usage=Rscript sampling_VOC_pg.R [metadata_file] [seqs_list_file] [country] [date_lower_limit] [date_upper_limit] [type_of_sampling] [VOC] [sample_size]")
  break
}
if (!require(dplyr, quietly = T)) install.packages("dplyr")
library(dplyr)
if (!require(tidyr, quietly = T)) install.packages("tidyr")
library(tidyr)
metadata<-read.csv(file=args[1],sep="\t",header=TRUE)
nexlist<-readLines(args[2])
f_d_c<-metadata%>%
  dplyr::filter(Country==args[3])%>%
  mutate(date=as.Date(date))%>%
  filter(date>args[4])%>%
  filter(date<args[5])%>%
  filter(name_analysis%in%nexlist)
#samples one by date####
p_s_1d<-function(voi,N_size){
  f_d_c_l<-f_d_c%>%
    dplyr::filter(grepl(voi,variant))%>%
    drop_na(date)
  f_d_c_l$date<-as.Date(f_d_c_l$date)
  f_d_c_l$City[f_d_c_l$City=="Amazona"]<-"Amazonas"
  f_d_c_l$City[f_d_c_l$City=="HuancavelIca"]<-"Huancavelica"
  f_d_c_l$City[f_d_c_l$City=="Huánuco"]<-"Huanuco"
  f_d_c_l$City[f_d_c_l$City=="Lalibertad"]<-"LaLibertad"
  f_d_c_l$City[f_d_c_l$City=="MadreDeDios"]<-"MadredeDios"
  f_d_c_l$region<-f_d_c_l$City
  f_d_c_l$region[f_d_c_l$City%in%c("Tumbes","Piura","Lambayeque","Cajamarca","Ancash","LaLibertad")]<-"North"
  f_d_c_l$region[f_d_c_l$City%in%c("Loreto","Amazonas","SanMartin")]<-"Northeast"
  f_d_c_l$region[f_d_c_l$City%in%c("Huanuco","Ucayali","Pasco","Junin")]<-"Mideast"
  f_d_c_l$region[f_d_c_l$City%in%c("Lima","Callao")]<-"Center"
  f_d_c_l$region[f_d_c_l$City%in%c("Ica","Arequipa","Moquegua","Tacna")]<-"South"
  f_d_c_l$region[f_d_c_l$City%in%c("Huancavelica","Ayacucho","Apurimac","Cusco","Puno","MadredeDios")]<-"Southeast"
  f_gr<-f_d_c_l%>%
    group_by(date,region)%>%
    summarise(count=n())
  var_list<-list()
  for (i in 1:nrow(f_gr)){
    row<-f_gr[i,] #crea obj de 1 fila
    var_list[[paste(row$date,row$region,"_")]]<-f_d_c_l%>%
      dplyr::filter(date==row$date,region==row$region)%>%
      dplyr::sample_n(size=1)
  }
  var<-do.call(rbind,var_list)
  if(nrow(var)>200){
    var<-var%>%
      sample_n(size=as.numeric(N_size))
  }
  names1<-paste(voi,"_s1_pg.list",sep="")
  tabname<-paste(voi,"_s1_pg.tsv",sep="")
  locname<-paste(voi,"_s1_region_pg.tsv",sep="")
  var_loc<-var%>%
    dplyr::select(name_analysis,region)
  writeLines(var$name_analysis,names1)
  write.table(var,tabname,sep="\t",quote=FALSE,row.names=FALSE)
  write.table(var_loc,locname,sep="\t",quote=FALSE,row.names=FALSE,col.names=F)
}
#random samples####
p_s_Nd<-function(voi,N_size){
  f_d_c_l<-f_d_c%>%
    dplyr::filter(grepl(voi,variant))%>%
    drop_na(date)
  f_d_c_l$date<-as.Date(f_d_c_l$date)
  f_d_c_l$City[f_d_c_l$City=="Amazona"]<-"Amazonas"
  f_d_c_l$City[f_d_c_l$City=="HuancavelIca"]<-"Huancavelica"
  f_d_c_l$City[f_d_c_l$City=="Huánuco"]<-"Huanuco"
  f_d_c_l$City[f_d_c_l$City=="Lalibertad"]<-"LaLibertad"
  f_d_c_l$City[f_d_c_l$City=="MadreDeDios"]<-"MadredeDios"
  f_d_c_l$region<-f_d_c_l$City
  f_d_c_l$region[f_d_c_l$City%in%c("Tumbes","Piura","Lambayeque","Cajamarca","Ancash","LaLibertad")]<-"North"
  f_d_c_l$region[f_d_c_l$City%in%c("Loreto","Amazonas","SanMartin")]<-"Northeast"
  f_d_c_l$region[f_d_c_l$City%in%c("Huanuco","Ucayali","Pasco","Junin")]<-"Mideast"
  f_d_c_l$region[f_d_c_l$City%in%c("Lima","Callao")]<-"Center"
  f_d_c_l$region[f_d_c_l$City%in%c("Ica","Arequipa","Moquegua","Tacna")]<-"South"
  f_d_c_l$region[f_d_c_l$City%in%c("Huancavelica","Ayacucho","Apurimac","Cusco","Puno","MadredeDios")]<-"Southeast"
  f_gr<-f_d_c_l%>%
    group_by(region)%>%
    summarise(count=n())%>%
    mutate(percent=count/sum(count))%>%
    mutate(sample=round(percent*as.numeric(N_size),0))
  var_list<-list()
  for (i in 1:nrow(f_gr)){
    row<-f_gr[i,] #crea obj de 1 fila
    var_list[[row$region]]<-f_d_c_l%>%
      dplyr::filter(region==row$region)%>%
      dplyr::sample_n(size=row$sample)
  }
  var<-do.call(rbind,var_list)
  if(nrow(var)>200){
    var<-var%>%
      sample_n(size=as.numeric(N_size))
  }
  names1<-paste(voi,"_sN_pg.list",sep="")
  tabname<-paste(voi,"_sN_pg.tsv",sep="")
  locname<-paste(voi,"_sN_region_pg.tsv",sep="")
  var_loc<-var%>%
    dplyr::select(name_analysis,region)
  writeLines(var$name_analysis,names1)
  write.table(var,tabname,sep="\t",quote=FALSE,row.names=FALSE)
  write.table(var_loc,locname,sep="\t",quote=FALSE,row.names=FALSE,col.names=F)
}
#processing####
if(args[6]=="s1"){
  p_s_1d(args[7],args[8]) 
}else{
  p_s_Nd(args[7],args[8]) 
}
