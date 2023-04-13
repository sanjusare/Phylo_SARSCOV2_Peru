#1:set working directory####
args<-commandArgs(trailingOnly = TRUE)
#2:load libraries and establish settings####
library(dplyr)
library(reshape2)
library(ggplot2)
library(openxlsx)
#3:read files####
d1<-read.table(args[1],sep="\t",header=T)
d2<-read.table(args[2],sep="\t",header=T)
d3<-read.table(args[3],sep="\t",header=T)
#4:Search of characteristic mutations####
d_f1<-d1%>%
  dplyr::select(Position,A,C,G,T,X.)%>%
  dplyr::mutate(Position=Position+203)%>%
  melt(id.vars="Position")%>%
  dplyr::mutate(code=paste(variable,Position,sep=""),
                sublineage="SubG1")
d_f2<-d2%>%
  dplyr::select(Position,A,C,G,T,X.)%>%
  dplyr::mutate(Position=Position+203)%>%
  melt(id.vars="Position")%>%
  dplyr::mutate(code=paste(variable,Position,sep=""),
                sublineage="SubG2")
d_f3<-d3%>%
  dplyr::select(Position,A,C,G,T,X.)%>%
  dplyr::mutate(Position=Position+203)%>%
  melt(id.vars="Position")%>%
  dplyr::mutate(code=paste(variable,Position,sep=""),
                sublineage="SubG3")
d_f1_ch<-d_f1%>%
  dplyr::filter(value>=0.8)
d_f2_ch<-d_f2%>%
  dplyr::filter(value>=0.8)
d_f3_ch<-d_f3%>%
  dplyr::filter(value>=0.8)

diff_1<-c()
for(i in 1:nrow(d_f1_ch)){
  if(d_f2$value[d_f2$code==d_f1_ch$code[i]]<=0.2&
     d_f3$value[d_f3$code==d_f1_ch$code[i]]<=0.2){
    diff_1<-append(diff_1,d_f1_ch$code[i])
  }
}

diff_2<-c()
for(i in 1:nrow(d_f2_ch)){
  if(#d_f1$value[d_f1$code==d_f2_ch$code[i]]<=0.2&
     d_f3$value[d_f3$code==d_f2_ch$code[i]]<=0.2){
    diff_2<-append(diff_2,d_f2_ch$code[i])
  }
}

diff_3<-c()
for(i in 1:nrow(d_f3_ch)){
  if(#d_f1$value[d_f1$code==d_f3_ch$code[i]]<=0.2&
     d_f3$value[d_f2$code==d_f3_ch$code[i]]<=0.2){
    diff_3<-append(diff_3,d_f3_ch$code[i])
  }
}

differences<-c(diff_1,diff_2,diff_3)
#5:Generating table with characteristic mutations of Lambda sublineages####
d_sel1<-d_f1_ch%>%
  filter(Position%in%as.numeric(gsub("[A-Z]","",differences)))
d_sel2<-d_f2_ch%>%
  filter(Position%in%as.numeric(gsub("[A-Z]","",differences)))
d_sel3<-d_f3_ch%>%
  filter(Position%in%as.numeric(gsub("[A-Z]","",differences)))
d_sel<-c(d_sel1$code,d_sel2$code,d_sel3$code)

d_self1<-d_f1%>%
  filter(code%in%d_sel)
d_self2<-d_f2%>%
  filter(code%in%d_sel)
d_self3<-d_f3%>%
  filter(code%in%d_sel)

d_f_p<-rbind(d_self1,d_self2,d_self3)%>%
  dplyr::select(Code=code,Relative_frequency=value,Sublineage=sublineage)%>%
  dplyr::arrange(Code)

write.xlsx(d_f_p,"Characteristic_subGamma_mutations.xlsx")
writeLines(unique(as.character(as.numeric(gsub("[A-Z]","",d_sel))-202)),"selected_Gamma_position.list")
