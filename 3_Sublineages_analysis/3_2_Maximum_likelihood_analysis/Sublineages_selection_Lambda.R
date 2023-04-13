#1:set working directory####
args<-commandArgs(trailingOnly = TRUE)
#2:load libraries and establish settings####
library(treeio)
library(ggtree)
library(dplyr)
library(ggplot2)
library(ggnewscale)
library(ggpubr)
theme_set(theme_bw(base_size=12))
Sys.setlocale("LC_TIME", "English")
colreg<-c("center"="firebrick","north"="dodgerblue","south"="forestgreen","mid-east"="darkmagenta","north-east"="goldenrod","south-east"="hotpink3")
collin<-c("SubL1"="#58578a","SubL2"="#8ad3db")
#3:read files####
#read the maximum likelihoood tree
t<-read.iqtree(args[1])
#read the metadata of genomes from Peru and filter those that are present in the tree
d<-read.csv(args[2],sep="\t")%>%
  dplyr::filter(name_analysis%in%t@phylo[["tip.label"]])%>%
  mutate(date=as.Date(date))
#read the file of the grouping information of cities in regions
reg<-read.csv(args[3],sep="\t")
#4:Sublineages identification####
#metadata arrangements
d$City[d$City=="Amazona"]<-"Amazonas"
d$City[d$City=="HuancavelIca"]<-"Huancavelica"
d$City[d$City=="Huánuco"]<-"Huanuco"
d$City[d$City=="Lalibertad"]<-"LaLibertad"
d$City[d$City=="MadreDeDios"]<-"MadredeDios"
d<-dplyr::left_join(d,reg)
#selection of sublineages
p_t<-ggtree(t,color="lightgrey")%<+%d
r_node<-p_t[["data"]]$parent[p_t[["data"]]$node==p_t[["data"]]$parent]
nodes<-p_t[["data"]]%>%
  dplyr::filter(isTip==FALSE&node!=r_node&SH_aLRT>=70&UFboot>=70)
sel_nodes<-c()
tips<-list()
for(i in 1:nrow(nodes)){
  print(i)
  t1<-tree_subset(t,node=nodes$node[i],levels_back=0)
  if(length(t1@phylo[["tip.label"]])>=400){
    tips[[as.character(nodes$node[i])]]<-t1@phylo[["tip.label"]]
    sel_nodes<-append(sel_nodes,nodes$node[i])
  }
}
sel_lins<-list()
for(i in 1:length(sel_nodes)){
  print(i)
  if(i!=length(sel_nodes)){
    test_tips<-tips[[i]][!(tips[[i]]%in%tips[[i+1]])]
  }else{
    test_tips<-tips[[i]]
  }
  if(length(test_tips)>=400){
    sel_lins[[as.character(sel_nodes[i])]]<-test_tips
  }
}
#save the tips that belong to each sublineage
for(i in 1:length(sel_lins)){
  name_lin<-paste("sublineage_SubL",i,".list",sep="")
  writeLines(sel_lins[[i]],name_lin)
}
#5:plotting the sublineages in the tree####
#preparing data
t_grouped<-t%>%
  groupOTU(sel_lins,group_name="SubLsT")
p_leg<-ggtree(t_grouped,aes(color=SubLsT),show.legend=F)%<+%d+
  scale_color_manual(values=c("lightgrey","#58578a","#8ad3db"))+
  new_scale_color()+
  geom_tippoint(aes(subset=label%in%unlist(sel_lins),color=Region))+
  scale_color_manual(values=colreg)+
  geom_nodelab(aes(subset=node==as.numeric(names(sel_lins)[1])),color="#58578a",nudge_x=-0.00009,nudge_y=20)+
  geom_nodepoint(aes(subset=node==as.numeric(names(sel_lins)[1])),color="#58578a",shape=15,size=4)+
  geom_nodelab(aes(subset=node==as.numeric(names(sel_lins)[2])),color="#8ad3db",nudge_x=-0.00009,nudge_y=20)+
  geom_nodepoint(aes(subset=node==as.numeric(names(sel_lins)[2])),color="#8ad3db",shape=15,size=4)+
  theme(legend.position="bottom")+
  guides(color=guide_legend(nrow=1))
leg<-get_legend(p_leg)
#specify that the data is from Lambda
t_grouped_lam<-t_grouped
d_lam<-d
sel_lins_lam<-sel_lins
#plotting the tree
p_t_fin_lam<-p_leg+
  theme(legend.position="none")
ggsave("Maximum_likelihood_Lambda.jpg",plot=p_t_fin_lam,dpi=500,height=5,width=5)



