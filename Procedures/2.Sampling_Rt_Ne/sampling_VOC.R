args<-commandArgs(trailingOnly = TRUE)
Sys.setlocale("LC_TIME", "English")
if(length(args)>7| length(args)<8){
  print("usage=Rscript sampling_VOC.R [metadata_file] [seqs_list_file] [country] [date_lower_limit] [date_upper_limit] [type_of_sampling] [VOC] [sample_size]")
  break
}
if (!require(dplyr, quietly = T)) install.packages("dplyr")
library(dplyr)
if (!require(tidyr, quietly = T)) install.packages("tidyr")
library(tidyr)
cols<-c("Omicron"="#db8a8a","Gamma"="#dbd48a","Lambda"="#8fdb8a","Delta"="#8ad3db","Mu"="#a18adb","Alpha"="#db8ac0","darkgrey")
metadata<-read.csv(file=args[1],sep="\t",header=TRUE)
nexlist<-readLines(args[2])
f_d_c<-metadata%>%
  dplyr::filter(Country==args[3])%>%
  mutate(date=as.Date(date))%>%
  filter(date>args[4])%>%
  filter(date<args[5])%>%
  filter(name_analysis%in%nexlist)
#samples one by date####
p_s_1d<-function(voi){
  f_d_c_l<-f_d_c%>%
    dplyr::filter(grepl(voi,variant))%>%
    drop_na(date)
  f_d_c_l$date<-as.Date(f_d_c_l$date)
  f_gr<-f_d_c_l%>%
    group_by(date)%>%
    summarise(count=n())
  var_list<-list()
  for (i in 1:nrow(f_gr)){
    row<-f_gr[i,] #crea obj de 1 fila
    var_list[[row$date]]<-f_d_c_l%>%
      dplyr::filter(date==row$date)%>%
      dplyr::sample_n(size=1)
  }
  var<-do.call(rbind,var_list)
  names1<-paste(voi,"_s1.list",sep="")
  tabname<-paste(voi,"_s1.tsv",sep="")
  writeLines(var$name_analysis,names1)
  write.table(var,tabname,sep="\t",quote=FALSE,row.names=FALSE)
  var_gr<-var%>%
    group_by(date)%>%
    summarise(count=n())
}
#random samples####
p_s_Nd<-function(voi,N_size){
  f_d_c_l<-f_d_c%>%
    dplyr::filter(grepl(voi,variant))%>%
    drop_na(date)
  f_d_c_l$date<-as.Date(f_d_c_l$date)
  f_gr<-f_d_c_l%>%
    dplyr::sample_n(size=N_size)
  names1<-paste(voi,"_sN.list",sep="")
  tabname<-paste(voi,"_sN.tsv",sep="")
  writeLines(f_gr$name_analysis,names1)
  write.table(f_gr,tabname,sep="\t",quote=FALSE,row.names=FALSE)
}
#processing####
if(args[6]=="s1"){
  p_s_1d(args[7]) 
}else{
  p_s_Nd(args[7],args[8]) 
}
