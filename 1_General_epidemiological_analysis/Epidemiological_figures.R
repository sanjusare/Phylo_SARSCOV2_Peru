#set working directory####
args<-commandArgs(trailingOnly = TRUE)
#load libraries and establish settings####
library(dplyr)
library(lubridate)
library(ggplot2)
library(tidyr)
library(zoo)
library(ggpubr)
library(ggtext)
theme_set(theme_bw(base_size=15))
Sys.setlocale("LC_TIME", "English")
cols<-c("Omicron"="#db8a8a","Gamma"="#dbd48a","Lambda"="#8fdb8a","Delta"="#8ad3db","Mu"="#a18adb","Alpha"="#db8ac0","Other"="darkgrey")
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
  labs(y="Reported cases")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(angle=90,hjust=0,vjust=0.5))+
  scale_x_datetime(date_breaks="1 month",date_labels="%Y-%b",limits=c(as.POSIXct("2020-01-30"),as.POSIXct("2022-04-30")))+
  scale_y_continuous(limits=c(0,45000),
                     breaks=seq(0,45000,4000),
                     sec.axis=sec_axis(trans=~./10,name="<span style='color:red'>Number of deaths</span>",breaks=seq(0,4500,400)))+
  theme(legend.position="none",
        axis.title.y=element_markdown(),
        axis.title.y.right=element_markdown(),
        axis.text.y.right=element_text(color="red"),
        axis.ticks.y.right=element_line(color="red"))+
  annotate("label",x=c(as.POSIXct("2020-07-15"),as.POSIXct("2021-04-15"),as.POSIXct("2022-01-25")),y=c(45000,45000,45000),label=c("1st Wave","2nd Wave","3rd wave"),fill="darkgrey",alpha=0.25)
ggsave("Reported_cases_deaths.jpg",plot=p_cas_dea,dpi=300,height=5,width=7.5)
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
                     sec.axis=sec_axis(trans=~./10,name=expression(Total~vaccinations~(log[10])),breaks=seq(0,10,1)))+
  theme(legend.position="none",
        axis.title.y.right=element_text(color="dodgerblue"),
        axis.text.y.right=element_text(color="dodgerblue"),
        axis.ticks.y.right=element_line(color="dodgerblue"))
ggsave("Stringency_index_Vaccinations.jpg",plot=p_SIVac,dpi=300,height=5,width=7.5)
#3IND_Figure_of_relative_prevalence_by_VOCI####
df_p<-tabla%>%
  dplyr::filter(Country=="Peru")%>%
  drop_na(date)
df_p$month<-as.numeric(sapply(strsplit(df_p$date,"-"),"[",2))
df_p$M<-as.character(month(ymd(010101) + months(df_p$month-1),label=TRUE,abbr=TRUE))
df_p$year<-sapply(strsplit(df_p$date,"-"),"[",1)
df_p$y_m<-paste(df_p$year,df_p$M,sep="-")
df_p$VOC<-sapply(strsplit(df_p$variant," "),"[",2)
df_p$VOC[is.na(df_p$VOC)]<-"Other"
for(i in 1:nrow(df_p)){
  if(!(df_p$VOC[i] %in% c("Lambda","Gamma","Alpha","Delta","Mu","Omicron","Other"))){
    df_p$VOC[i]<-"Other"
  }
}
df_p_cg<-df_p%>%dplyr::group_by(y_m,VOC)%>%
  dplyr::summarise(count=n())%>%
  dplyr::mutate(percentage=count/sum(count))
df_p_cg$y_m<-factor(df_p_cg$y_m,levels=c("2020-Feb","2020-Mar","2020-Apr","2020-May","2020-Jun","2020-Jul","2020-Aug","2020-Sep","2020-Oct","2020-Nov","2020-Dec","2021-Jan","2021-Feb","2021-Mar","2021-Apr","2021-May","2021-Jun","2021-Jul","2021-Aug","2021-Sep","2021-Oct","2021-Nov","2021-Dec","2022-Jan","2022-Feb","2022-Mar","2022-Apr"))
df_p_cg$y_m<-as.Date(paste("1-",df_p_cg$y_m,sep=""),"%d-%Y-%B")
df_p_cg$VOC<-factor(df_p_cg$VOC,levels=c("Lambda","Gamma","Alpha","Delta","Mu","Omicron","Other"))
p_rel_prev<-ggplot(subset(df_p_cg),aes(x=as.POSIXct(y_m),y=percentage,fill=VOC))+
  geom_bar(stat="identity")+
  theme(axis.text.x=element_text(angle=90,hjust=0,vjust=0.5),axis.title.x=element_blank())+
  labs(y="Relative prevalence")+
  scale_fill_manual(values=cols)+
  theme(legend.title=element_blank(),
        legend.position="bottom")+
  guides(fill=guide_legend(nrow=1))+
  scale_x_datetime(date_breaks="1 month",date_labels="%Y-%b",limits=c(as.POSIXct("2020-01-30"),as.POSIXct("2022-04-30")))
ggsave("Relative_Prevalence_VOCI.jpg",plot=p_rel_prev,dpi=500,height=5,width=7.5)
#4depon1,2,3_Integrated_Figure####
leg<-get_legend(p_rel_prev)
p_rel_prev_woleg<-p_rel_prev+
  theme(legend.position="none")
Figure1<-ggarrange(p_cas_dea,p_SIVac,p_rel_prev_woleg,ncol=1,nrow=3,align="hv",labels=c("A","B","C"),common.legend=T,legend.grob=leg,legend="bottom")
ggsave(filename="General_analysis.jpg",plot=Figure1,dpi=300,height=11,width=7.5)
#5IND_Figure_of_tests_per_inhabitant_and_deaths_by_country####
#arrange data
dea_tests<-cas_wor%>%
  dplyr::filter(location%in%c("Peru","Argentina","Chile"))%>%
  dplyr::select(location,Date=date,new_deaths,new_cases,new_tests_smoothed_per_thousand)%>%
  dplyr::mutate(new_deaths=replace_na(new_deaths,0),
                new_cases=replace_na(new_cases,0),
                new_tests_per_thousand=replace_na(new_tests_smoothed_per_thousand,0))%>%
  dplyr::mutate(Date=as.Date(Date))%>%
  dplyr::mutate(rolldeath_14=zoo::rollmean(new_deaths,k=14,fill=0),
                rollcases_14=zoo::rollmean(new_cases,k=14,fill=0))
#plot deaths and cases by country
p_dea<-ggplot(data=dea_tests[7:(nrow(dea_tests)-7),],aes(x=as.POSIXct(Date)))+
  geom_rect(aes(xmin=as.POSIXct("2021-01-01"),xmax=as.POSIXct("2021-08-01"),ymin=0,ymax=32000),fill="gray85",alpha=0.2)+
  geom_line(aes(y=round(rolldeath_14,0)*20,color="red"))+
  geom_line(aes(y=round(rollcases_14,0)))+
  labs(y="Reported cases")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(angle=90,hjust=0,vjust=0.5))+
  scale_x_datetime(date_breaks="1 month",date_labels="%Y-%b",limits=c(as.POSIXct("2020-12-01"),as.POSIXct("2021-09-01")))+
  scale_y_continuous(limits=c(0,33000),
                     breaks=seq(0,33000,3000),
                     sec.axis=sec_axis(trans=~./20,name="<span style='color:red'>Number of deaths</span>",breaks=seq(0,1650,165)))+
  facet_wrap(~location,ncol=1,nrow=3)+
  theme(legend.position="none",
        axis.title.y=element_markdown(),
        axis.title.y.right=element_markdown(),
        axis.text.y.right=element_text(color="red"),
        axis.ticks.y.right=element_line(color="red"))+
  annotate("label",x=c(as.POSIXct("2021-04-15")),y=c(33000),label=c("2nd Wave"),fill="darkgrey",alpha=0.25)
ggsave("Reported_deaths_cases_country.jpg",plot=p_dea,dpi=300,height=10,width=5)
#plot tests per thousand inhabitants
p_tests<-ggplot(data=dea_tests[7:(nrow(dea_tests)-7),],aes(x=as.POSIXct(Date)))+
  geom_rect(aes(xmin=as.POSIXct("2021-01-01"),xmax=as.POSIXct("2021-08-01"),ymin=0,ymax=5),fill="gray85",alpha=0.2)+
  geom_line(aes(y=new_tests_per_thousand))+
  labs(y="Tests per thousand")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(angle=90,hjust=0,vjust=0.5))+
  scale_x_datetime(date_breaks="1 month",date_labels="%Y-%b",limits=c(as.POSIXct("2020-12-01"),as.POSIXct("2021-09-01")))+
  scale_y_continuous(position="right")+
  facet_wrap(~location,ncol=1,nrow=3)+
  theme(legend.position="none")+
  annotate("label",x=c(as.POSIXct("2021-04-15")),y=c(6),label=c("2nd Wave"),fill="darkgrey",alpha=0.25)
ggsave("Test_by_inhabitant.jpg",plot=p_tests,dpi=300,height=10,width=5)
#plotting together deaths,cases and tests
p_dea_tests<-ggarrange(p_dea,p_tests,labels=c("A","B"),ncol=2,nrow=1,align="hv")
ggsave("Reported_deaths_cases_tests_country.jpg",plot=p_dea_tests,dpi=500,height=10,width=10)
