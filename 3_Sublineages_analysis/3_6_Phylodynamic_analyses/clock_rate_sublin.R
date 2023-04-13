#1:set working directory####
args<-commandArgs(trailingOnly = TRUE)
#2:load libraries and establish settings####
library(ggplot2)
library(dplyr)
theme_set(theme_bw(base_size=15))
Sys.setlocale("LC_TIME", "English")
#3:Read files####
files<-dir(pattern="*f2.log")
t_list<-list()
for(i in 1:length(files)){
  t<-read.csv(files[i],sep="\t",comment.char="#")
  lin<-sapply(strsplit(files[i],"_"),"[",1)
  t_list[[paste(lin,samp,sep="_")]]<-t%>%
    dplyr::select(clock_rate=names(t)[grepl("clockRate",names(t))])%>%
    dplyr::mutate(lineage=lin)
}
t<-do.call(rbind,t_list)
t$lineage[t$lineage=="SubL"]<-"SubL1+L2"
#4:Plot of clock rate distribution####
p_clock<-ggplot(t,aes(x=clock_rate))+
  geom_density(alpha=0.3)+
  labs(x="Clock rate (s/s/y)")+
  theme(axis.title.y=element_blank(),
        legend.position="bottom")+
  facet_wrap(~lineage)
ggsave("Clock_rate_sublineages.jpg",plot=p_clock,dpi=500,height=5,width=5)  
