#1:set working directory####
args<-commandArgs(trailingOnly = TRUE)
#2:load libraries and establish settings####
library(ggtree)
library(treeio)
library(ggplot2)
library(ggpubr)
theme_set(theme_bw(base_size=15))
Sys.setlocale("LC_TIME", "English")
colreg<-c("center"="firebrick","north"="dodgerblue","south"="forestgreen","mid-east"="darkmagenta","north-east"="goldenrod","south-east"="hotpink3")
collin<-c("SubL1+L2"="forestgreen","SubG1"="dodgerblue","SubG2"="firebrick","SubG3"="#ad71ab")
#3:read files####
files<-dir(pattern="*_MCT.tree")
#4:Extract the probabilities of each location to be the origin and the date of origin####
height_list<-list()
prob_list<-list()
for(i in 1:length(files)){
  tr<-read.beast(files[i])
  data<-data.frame(tip=tr@phylo[["tip.label"]],
                   date=as.Date(sapply(strsplit(tr@phylo[["tip.label"]],"\\|"),"[",4)))
  mrt<-Date2decimal(max(data$date))
  locset<-tr@data$location.set[tr@data$node==length(tr@phylo[["tip.label"]])+1][[1]]
  locprob<-tr@data$location.set.prob[tr@data$node==length(tr@phylo[["tip.label"]])+1][[1]]
  height_medi<-as.numeric(tr@data$CAheight_median[tr@data$node==length(tr@phylo[["tip.label"]])+1][[1]])
  height_mean<-as.numeric(tr@data$CAheight_mean[tr@data$node==length(tr@phylo[["tip.label"]])+1][[1]])
  height_low95<-as.numeric(tr@data$CAheight_0.95_HPD[tr@data$node==length(tr@phylo[["tip.label"]])+1][[1]][1])
  height_hig95<-as.numeric(tr@data$CAheight_0.95_HPD[tr@data$node==length(tr@phylo[["tip.label"]])+1][[1]][2])
  date_medi<-decimal2Date(mrt-height_medi)
  date_mean<-decimal2Date(mrt-height_mean)
  date_low95<-decimal2Date(mrt-height_low95)
  date_hig95<-decimal2Date(mrt-height_hig95)
  height_list[[as.character(i)]]<-data.frame(height_median=height_medi,
                                             height_mean=height_mean,
                                             height_low95=height_low95,
                                             height_hig95=height_hig95,
                                             date_median=date_medi,
                                             date_mean=date_mean,
                                             date_low95=date_low95,
                                             date_hig95=date_hig95,
                                             sample=as.character(sapply(strsplit(files[i],"_"),"[",1)),
                                             subsample=as.numeric(gsub("v2","",sapply(strsplit(files[i],"_"),"[",4))),
                                             mrt=mrt)
  prob_list[[as.character(i)]]<-data.frame(country=locset,
                                           prob=locprob,
                                           sample=as.character(sapply(strsplit(files[i],"_"),"[",1)),
                                           subsample=as.numeric(gsub("v2","",sapply(strsplit(files[i],"_"),"[",4))))
}
heights<-do.call(rbind,height_list)
probs<-do.call(rbind,prob_list)
probs$sample[probs$sample=="SubL"]<-"SubL1+L2"
heights$sample[heights$sample=="SubL"]<-"SubL1+L2"
#5:Calculate means by repetitions of location origin probabilities####
over_prob<-probs%>%
  dplyr::group_by(country,sample)%>%
  dplyr::summarise(prob=mean(prob))
over_prob$sample<-factor(over_prob$sample,levels=c("SubL1+L2","SubG1","SubG2","SubG3"),ordered=T)
#6:Calculate means by repetitions of origin dates probabilities####
over_heights<-heights%>%
  dplyr::group_by(sample)%>%
  dplyr::summarise(height_median=mean(height_median),
                   height_mean=mean(height_mean),
                   height_low95=mean(height_low95),
                   height_hig95=mean(height_hig95),
                   date_median=mean(date_median),
                   date_mean=mean(date_mean),
                   date_low95=mean(date_low95),
                   date_hig95=mean(date_hig95))
over_heights$sample<-factor(over_heights$sample,levels=c("SubL1+L2","SubG1","SubG2","SubG3"),ordered=T)
#7:Plots of probability of each region to be the origin of sublineages and date distribution of the origin####
pprobs<-ggplot(over_prob,aes(x=as.factor(sample),y=prob,fill=country))+
  geom_bar(stat="identity",position="dodge")+
  labs(y="Location probability of root")+
  scale_fill_manual(values=colreg)+
  theme(legend.title=element_blank(),
        legend.position=c(0.6,0.85),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=10),
        legend.key.size=unit(0.1,"in"),
        legend.margin=margin(0.01,0.01,0.01,0.01),
        legend.box.margin=margin(0.001,0.001,0.001,0.001),
        legend.box.spacing=unit(0.1,"in"),
        legend.text=element_text(size=7))+
  ylim(c(0,1))+
  guides(fill=guide_legend(ncol=1))
pheights<-ggplot(over_heights,aes(x=as.factor(sample),y=as.POSIXct(date_median)))+
  geom_pointrange(aes(ymin=as.POSIXct(date_low95),ymax=as.POSIXct(date_hig95),color=sample))+
  geom_point(aes(y=as.POSIXct(date_mean),fill=sample),shape=21,size=5,alpha=0.5)+
  labs(y="Root Date")+
  scale_y_datetime(date_breaks="1 month",date_labels="%Y-%b",limits=c(as.POSIXct("2020-04-01"),as.POSIXct("2021-06-01")))+
  theme(axis.title.x=element_blank(),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=10),
        legend.position="none")+
  scale_color_manual(values=collin)+
  scale_fill_manual(values=collin)
p_origin_lin<-ggarrange(pprobs,pheights,ncol=2,nrow=1,labels=c("A","B"),align="h",widths=c(0.9,1))
ggsave("Date_Location_of_Sublineages_origin.jpg",plot=p_origin_lin,dpi=500,height=4,width=8)
