#1:set working directory####
args<-commandArgs(trailingOnly = TRUE)
#2:load libraries and establish settings####
library(ggplot2)
library(ggtree)
library(treeio)
library(ggpubr)
theme_set(theme_bw(base_size=15))
Sys.setlocale("LC_TIME", "English")
colors<-c("#cbcbf7","#f7cbcb","#f7f4cb","thistle1","#91e3db","#abc4d9","#9c9181","aquamarine4","#a9c7af","#99a5d6","#f4b771","#ec74dc","#9ed5ff","#dedede")
count_colors<-c("Argentina"="dodgerblue","Chile"="forestgreen","Peru"="firebrick","Colombia"="orangered","Ecuador"="goldenrod","Mexico"="darkmagenta")
#3:figure of the Maximum Clade Credibility tree####
#arrange the data
data<-read.csv(args[3],sep="\t")
data$date<-as.Date(data$date)
data_order<-data%>%
  arrange(date)
#read MCC trees 
tar<-read.beast(args[1])
tpe<-read.beast(args[2])
#read summary of transitions
map_sum<-read.csv(args[4],sep="\t")
#generate the base image of the tree
p_tpe<-ggtree(tpe,mrsd=data_order$date[nrow(data_order)],color="lightgrey")%<+%data_order
p_tar<-ggtree(tar,mrsd=data_order$date[nrow(data_order)],color="lightgrey")%<+%data_order
#select the most ancient and high supported nodes from each country
pernodes_list<-list()
for(i in levels(as.factor(p_tpe[["data"]]$Country))){
  pernodes_list[[i]]<-p_tpe[["data"]]%>%
    dplyr::filter(location==i)%>%
    dplyr::mutate(CAheight_median=as.numeric(CAheight_median))%>%
    dplyr::filter(CAheight_median==max(CAheight_median,na.rm=T))
}
pernodes<-do.call(rbind,pernodes_list)
argnodes_list<-list()
for(i in levels(as.factor(p_tar[["data"]]$Country))){
  argnodes_list[[i]]<-p_tar[["data"]]%>%
    dplyr::filter(location==i)%>%
    dplyr::mutate(CAheight_median=as.numeric(CAheight_median))%>%
    dplyr::filter(CAheight_median==max(CAheight_median,na.rm=T))
}
argnodes<-do.call(rbind,argnodes_list)
#plot the trees
#Peru
rootper<-Date2decimal(max(as.Date(p_tpe[["data"]]$date),na.rm=T))-as.numeric(p_tpe[["data"]]$height[p_tpe[["data"]]$node==p_tpe[["data"]]$parent])
receper<-Date2decimal(max(as.Date(p_tpe[["data"]]$date),na.rm=T))
p_tpe2<-p_tpe+
  theme_tree2()+  
  geom_nodepoint(aes(color=location),size=1)+
  geom_nodepoint(aes(subset=node%in%pernodes$node,color=location),size=5,shape=21,stroke=1.5)+
  scale_color_manual(values=count_colors[names(count_colors)%in%c("Argentina","Peru","Chile")])+
  ggnewscale::new_scale_color()+
  geom_tippoint(aes(color=location),size=1,shape=15)+
  scale_color_manual(values=count_colors[names(count_colors)%in%c("Argentina","Peru","Chile")])+
  theme(axis.text.x=element_text(angle=90,hjust=0,vjust=0,size=10),
        panel.grid.major=element_line(size=0.5),
        legend.position="none",
        legend.title=element_blank(),
        legend.box="vertical",
        legend.direction="vertical")+
  ylim(c(-0.001,(length(tpe@phylo[["tip.label"]])+0.001)))+
  scale_x_continuous(limits=c(rootper-0.08,receper),breaks=seq(2020,2022,0.1),labels=decimal2Date(seq(2020,2022,0.1)))+
  annotate(geom="label",x=pernodes$x[pernodes$location=="Chile"]+0.20,y=pernodes$y[pernodes$location=="Chile"],color="forestgreen",alpha=0.5,
           size=3,
           label=paste(map_sum$q05date2[map_sum$locroot=="Peru"&map_sum$trans=="Peru to Chile"],
                       "-",
                       map_sum$q95date2[map_sum$locroot=="Peru"&map_sum$trans=="Peru to Chile"]))+
  annotate(geom="label",x=pernodes$x[pernodes$location=="Argentina"]+0.20,y=pernodes$y[pernodes$location=="Argentina"],color="dodgerblue",alpha=0.5,
           size=3,
           label=paste(map_sum$q05date2[map_sum$locroot=="Peru"&map_sum$trans=="Peru to Argentina"],
                       "-",
                       map_sum$q95date2[map_sum$locroot=="Peru"&map_sum$trans=="Peru to Argentina"]))+
  geom_taxalink(taxa1=pernodes$node[pernodes$location=="Peru"],taxa2=pernodes$node[pernodes$location=="Chile"],arrow=arrow(length=unit(0.2,"inches")),color="forestgreen")+
  geom_taxalink(taxa1=pernodes$node[pernodes$location=="Peru"],taxa2=pernodes$node[pernodes$location=="Argentina"],arrow=arrow(length=unit(0.2,"inches")),color="dodgerblue")
#argentina
rootarg<-Date2decimal(max(as.Date(p_tar[["data"]]$date),na.rm=T))-as.numeric(p_tar[["data"]]$height[p_tar[["data"]]$node==p_tar[["data"]]$parent])
recearg<-Date2decimal(max(as.Date(p_tar[["data"]]$date),na.rm=T))
p_tar2<-p_tar+
  theme_tree2()+
  geom_nodepoint(aes(color=location),size=1)+
  geom_nodepoint(aes(subset=node%in%argnodes$node,color=location),size=5,shape=21,stroke=1.5)+
  scale_color_manual(values=count_colors[names(count_colors)%in%c("Argentina","Peru","Chile")])+
  ggnewscale::new_scale_color()+
  geom_tippoint(aes(color=location),size=1,shape=15)+
  scale_color_manual(values=count_colors[names(count_colors)%in%c("Argentina","Peru","Chile")])+
  theme(axis.text.x=element_text(angle=90,hjust=0,vjust=0,size=10),
        panel.grid.major=element_line(size=0.5),
        legend.position="none",
        legend.title=element_blank(),
        legend.box="vertical",
        legend.direction="vertical")+
  ylim(c(-0.001,(length(tpe@phylo[["tip.label"]])+0.001)))+
  scale_x_continuous(limits=c(rootarg-0.08,recearg),breaks=seq(2020,2022,0.1),labels=decimal2Date(seq(2020,2022,0.1)))+
  annotate(geom="label",x=argnodes$x[argnodes$location=="Chile"]+0.20,y=argnodes$y[argnodes$location=="Chile"],color="forestgreen",alpha=0.5,
           size=3,
           label=paste(map_sum$q05date2[map_sum$locroot=="Argentina"&map_sum$trans=="Peru to Chile"],
                       "-",
                       map_sum$q95date2[map_sum$locroot=="Argentina"&map_sum$trans=="Peru to Chile"]))+
  annotate(geom="label",x=argnodes$x[argnodes$location=="Peru"]+0.20,y=argnodes$y[argnodes$location=="Peru"],color="firebrick",alpha=0.5,
           size=3,
           label=paste(map_sum$q05date2[map_sum$locroot=="Argentina"&map_sum$trans=="Argentina to Peru"],
                       "-",
                       map_sum$q95date2[map_sum$locroot=="Argentina"&map_sum$trans=="Argentina to Peru"]))+
  geom_taxalink(taxa1=argnodes$node[argnodes$location=="Argentina"],taxa2=argnodes$node[argnodes$location=="Peru"],arrow=arrow(length=unit(0.2,"inches")),color="firebrick")+
  geom_taxalink(taxa1=argnodes$node[argnodes$location=="Peru"],taxa2=argnodes$node[argnodes$location=="Chile"],arrow=arrow(length=unit(0.2,"inches")),color="forestgreen")
#4:Production of a composed figure####
rootpe<-ggarrange(p_tpe2,ncol=1,nrow=1)%>%
  annotate_figure(left=text_grob("Peru",rot=90))
rootar<-ggarrange(p_tar2,ncol=1,nrow=1)%>%
  annotate_figure(left=text_grob("Argentina",rot=90))
pfin<-ggarrange(rootpe,rootar,ncol=1,nrow=2)%>%
  annotate_figure(left=text_grob("Location root",rot=90))
ggsave("MCC_trees_of_Lambda_origin.jpg",plot=pfin,dpi=500,height=7.5,width=7)
