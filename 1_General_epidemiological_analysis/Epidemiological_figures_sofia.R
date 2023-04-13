#set working directory####
args<-commandArgs(trailingOnly = TRUE)
args<-c("ordered_metadata.tsv","owid-covid-data.csv")
#load libraries and establish settings####
library(dplyr)
library(lubridate)
library(ggplot2)
library(tidyr)
library(zoo)
library(ggpubr)
library(ggtext)
theme_set(theme_bw(base_size=15))
Sys.setlocale("LC_TIME", "Spanish")
cols<-c("Omicron"="#db8a8a","Gamma"="#dbd48a","Lambda"="#8fdb8a","Delta"="#8ad3db","Mu"="#a18adb","Alpha"="#db8ac0","Otros"="darkgrey")
colors<-c("#cbcbf7","#f7cbcb","#f7f4cb","thistle1","#91e3db","#abc4d9","#9c9181","aquamarine4","#a9c7af","#99a5d6","#f4b771","#ec74dc","#9ed5ff","#dedede")
#Read input files####
#obtained from the output of Preprocessing_metadata.R script
#"ordered_metadata.tsv"
tabla<-read.csv(args[1],sep="\t")
#download from https://github.com/owid/covid-19-data/tree/master/public/data
#"owid-covid-data.csv"
cas_wor<-read.csv(args[2])
#1IND_Figure_of_reported_cases_and_deaths####
cas_dea<-cas_wor%>%
  dplyr::filter(location=="Peru")%>%
  dplyr::select(Date=date,new_deaths,new_cases)%>%
  dplyr::mutate(new_deaths=replace_na(new_deaths,0),
                new_cases=replace_na(new_cases,0),)%>%
  dplyr::mutate(Date=as.Date(Date))%>%
  dplyr::mutate(rolldeath_14=zoo::rollmean(new_deaths,k=14,fill=0),
                rollcases_14=zoo::rollmean(new_cases,k=14,fill=0))
p_cas_dea<-ggplot(data=cas_dea[7:(nrow(cas_dea)-7),],aes(x=as.POSIXct(Date)))+
  geom_rect(aes(xmin=as.POSIXct("2020-04-01"),xmax=as.POSIXct("2020-11-01"),ymin=0,ymax=44000),fill="gray85",alpha=0.2)+
  geom_rect(aes(xmin=as.POSIXct("2021-01-01"),xmax=as.POSIXct("2021-08-01"),ymin=0,ymax=44000),fill="gray85",alpha=0.2)+
  geom_rect(aes(xmin=as.POSIXct("2021-12-01"),xmax=as.POSIXct("2022-04-01"),ymin=0,ymax=44000),fill="gray85",alpha=0.2)+
  geom_line(aes(y=round(rolldeath_14,0)*10,color="red"))+
  geom_line(aes(y=round(rollcases_14,0)))+
  labs(y="Casos Reportados")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(angle=90,hjust=0,vjust=0.5))+
  scale_x_datetime(date_breaks="1 month",date_labels="%Y-%b",limits=c(as.POSIXct("2020-01-30"),as.POSIXct("2022-04-30")))+
  scale_y_continuous(limits=c(0,45000),
                     breaks=seq(0,45000,4000),
                     sec.axis=sec_axis(trans=~./10,name="<span style='color:red'>Numero de descesos</span>",breaks=seq(0,4500,400)))+
  theme(legend.position="none",
        axis.title.y=element_markdown(),
        axis.title.y.right=element_markdown(),
        axis.text.y.right=element_text(color="red"),
        axis.ticks.y.right=element_line(color="red"))+
  annotate("label",x=c(as.POSIXct("2020-07-15"),as.POSIXct("2021-04-15"),as.POSIXct("2022-01-25")),y=c(45000,45000,45000),label=c("1era ola","2da ola","3era ola"),fill="darkgrey",alpha=0.25)
ggsave("Reported_cases_deaths_sofia.jpg",plot=p_cas_dea,dpi=300,height=5,width=7.5)
#2Depon1_Figure_of_stringency_index_and_vaccinations####
SI_VA_tab<-cas_wor%>%
  dplyr::filter(location=="Peru")%>%
  dplyr::select(Date=date,
                SI=stringency_index,
                Vac=total_vaccinations)%>%
  dplyr::mutate(Date=as.Date(Date),
                log_vac=log10(Vac))
p_SIVac<-ggplot(data=cas_dea[7:(nrow(cas_dea)-7),],aes(x=as.POSIXct(Date)))+
  geom_bar(data=SI_VA_tab,aes(y=SI),stat="identity",fill="gray85")+
  geom_line(data=SI_VA_tab,aes(y=log_vac*10),linetype="longdash",color="dodgerblue",size=1)+
  labs(y="Stringency index")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(angle=90,hjust=0,vjust=0.5))+
  scale_x_datetime(date_breaks="1 month",date_labels="%Y-%b",limits=c(as.POSIXct("2020-01-30"),as.POSIXct("2022-04-30")))+
  scale_y_continuous(limits=c(0,100),
                     breaks=seq(0,100,10),
                     sec.axis=sec_axis(trans=~./10,name=expression(Numero~de~vacunaciones~(log[10])),breaks=seq(0,10,1)))+
  theme(legend.position="none",
        axis.title.y.right=element_text(color="dodgerblue"),
        axis.text.y.right=element_text(color="dodgerblue"),
        axis.ticks.y.right=element_line(color="dodgerblue"))
ggsave("Stringency_index_Vaccinations_sofia.jpg",plot=p_SIVac,dpi=300,height=5,width=7.5)
#3IND_Figure_of_relative_prevalence_by_VOCI####
Sys.setlocale("LC_TIME", "English")
df_p<-tabla%>%
  dplyr::filter(Country=="Peru")%>%
  drop_na(date)
df_p$month<-as.numeric(sapply(strsplit(df_p$date,"-"),"[",2))
df_p$M<-as.character(month(ymd(010101) + months(df_p$month-1),label=TRUE,abbr=TRUE))
df_p$year<-sapply(strsplit(df_p$date,"-"),"[",1)
df_p$y_m<-paste(df_p$year,df_p$M,sep="-")
df_p$VOC<-sapply(strsplit(df_p$variant," "),"[",2)
df_p$VOC[is.na(df_p$VOC)]<-"Otros"
for(i in 1:nrow(df_p)){
  if(!(df_p$VOC[i] %in% c("Lambda","Gamma","Alpha","Delta","Mu","Omicron","Otros"))){
    df_p$VOC[i]<-"Otros"
  }
}
df_p_cg<-df_p%>%dplyr::group_by(y_m,VOC)%>%
  dplyr::summarise(count=n())%>%
  dplyr::mutate(percentage=count/sum(count))
df_p_cg$y_m<-as.Date(paste("1-",df_p_cg$y_m,sep=""),"%d-%Y-%B")
df_p_cg$VOC<-factor(df_p_cg$VOC,levels=c("Lambda","Gamma","Alpha","Delta","Mu","Omicron","Otros"))
Sys.setlocale("LC_TIME", "Spanish")
p_rel_prev<-ggplot(subset(df_p_cg),aes(x=as.POSIXct(y_m),y=percentage,fill=VOC))+
  geom_bar(stat="identity")+
  theme(axis.text.x=element_text(angle=90,hjust=0,vjust=0.5),axis.title.x=element_blank())+
  labs(y="Prevalencia relativa")+
  scale_fill_manual(values=cols)+
  theme(legend.title=element_blank(),
        legend.position="bottom")+
  guides(fill=guide_legend(nrow=1))+
  scale_x_datetime(date_breaks="1 month",date_labels="%Y-%b",limits=c(as.POSIXct("2020-01-30"),as.POSIXct("2022-04-30")))
ggsave("Relative_Prevalence_VOCI_sofia.jpg",plot=p_rel_prev,dpi=500,height=5,width=7.5)
#4depon1,2,3_Integrated_Figure####
leg<-get_legend(p_rel_prev)
p_rel_prev_woleg<-p_rel_prev+
  theme(legend.position="none")
Figure1<-ggarrange(p_cas_dea,p_SIVac,p_rel_prev_woleg,ncol=1,nrow=3,align="hv",labels=c("A","B","C"),common.legend=T,legend.grob=leg,legend="bottom")
ggsave(filename="General_analysis.jpg",plot=Figure1,dpi=300,height=11,width=7.5)
