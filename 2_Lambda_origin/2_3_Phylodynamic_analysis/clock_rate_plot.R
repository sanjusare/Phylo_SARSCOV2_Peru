#1:set working directory####
args<-commandArgs(trailingOnly = TRUE)
#2:load libraries and establish settings####
library(ggplot2)
library(dplyr)
theme_set(theme_bw(base_size=15))
Sys.setlocale("LC_TIME", "English")
#3:Read files####
files<-args[1]
t_list<-list()
for(i in 1:length(files)){
  t<-read.csv(files[i],sep="\t",comment.char="#")
  t_list[[as.character(i)]]<-t%>%
    dplyr::select(clock_rate=names(t)[grepl("clockRate",names(t))])%>%
    dplyr::mutate(sample=i)
}
t<-do.call(rbind,t_list)
#4:Plot of clock rate distribution####
p_clock<-ggplot(t,aes(x=clock_rate))+
  geom_density(alpha=0.3,fill="forestgreen")+
  labs(x="Clock rate (s/s/y)")+
  theme(axis.title.y=element_blank(),
        legend.position="none")
ggsave("Clock_rate_distribution.jpg",plot=p_clock,dpi=500,height=5,width=5)  
