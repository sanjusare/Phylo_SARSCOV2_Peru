args<-commandArgs(trailingOnly = TRUE)
Sys.setlocale("LC_TIME", "English")
if(length(args)>4| length(args)<4){
  print("usage=Rscript select_new_lins_from_tree.R [metadata_file] [seqs_list_file] [country] [VOC]")
  break
}
if (!require(dplyr, quietly = T)) install.packages("dplyr")
library(dplyr)
if (!require(tidyr, quietly = T)) install.packages("tidyr")
library(tidyr)
if (!require(ggtree, quietly = T)) install.packages("ggtree")
library(ggtree)
if (!require(ggplot2, quietly = T)) install.packages("ggplot2")
library(ggplot2)
if (!require(ggpubr, quietly = T)) install.packages("ggpubr")
library(ggpubr)
if (!require(ggnewscale, quietly = T)) install.packages("ggnewscale")
library(ggnewscale)
#searching well supported peruvian clades####
data<-read.table(args[1],sep="\t",quote="",header=TRUE)
data$name_analysis<-gsub("'","_",data$name_analysis)
data$Peru<-"no"
data$Peru[data$Country==args[3]]<-"yes"
data$date<-as.Date(data$date)
t_f_iq<-read.iqtree(args[2])
t_f_iq@phylo[["tip.label"]]<-sapply(strsplit(t_f_iq@phylo[["tip.label"]],"\\|"),"[",2)
data_clean<-data%>%
  dplyr::filter(EPI_ISL%in%t_f_iq@phylo[["tip.label"]])%>%
  dplyr::select(!name_analysis)
p_t_f_iq<-ggtree(t_f_iq,color="lightgrey")%<+%data_clean
nodes<-p_t_f_iq[["data"]]%>%
  dplyr::filter(Country==args[3])
r_node<-p_t_f_iq[["data"]]$parent[p_t_f_iq[["data"]]$node==p_t_f_iq[["data"]]$parent]
saved_nodes_list<-list()
for(z in nodes$node){
  saved_nodes_list[[as.character(z)]]<-c(p_t_f_iq[["data"]]$parent[p_t_f_iq[["data"]]$node==z])
  for(i in 1:1000000){
    if(saved_nodes_list[[as.character(z)]][i]==r_node){
      break
    }else{
      saved_nodes_list[[as.character(z)]]<-append(saved_nodes_list[[as.character(z)]],p_t_f_iq[["data"]]$parent[p_t_f_iq[["data"]]$node==saved_nodes_list[[as.character(z)]][i]])
    }
  }
}
saved_nodes<-unlist(saved_nodes_list)
saved_nodes_unique<-unique(saved_nodes)
p_t_f_iq[["data"]]$SH_aLRT[is.na(p_t_f_iq[["data"]]$SH_aLRT)]<-0
p_t_f_iq[["data"]]$UFboot[is.na(p_t_f_iq[["data"]]$UFboot)]<-0
n_sel_post<-c()
for(i in saved_nodes_unique){
  if(i==r_node){
  }else if((p_t_f_iq[["data"]]$SH_aLRT[p_t_f_iq[["data"]]$node==i])>=90 &&
           (p_t_f_iq[["data"]]$UFboot[p_t_f_iq[["data"]]$node==i])>=90){
    n_sel_post<-append(n_sel_post,i)
  }
}
sel_trees<-list()
sel_trees_tips_lists<-list()
for(i in 1:length(n_sel_post)){
  t1<-tree_subset(t_f_iq,node=n_sel_post[i],levels_back=0)
  d_tip<-data_clean%>%
    dplyr::filter(EPI_ISL%in%t1@phylo[["tip.label"]])
  d_t_p<-d_tip%>%
    dplyr::filter(Country==args[3])
  if(nrow(d_t_p)>=round(nrow(nodes)*0.01,0)&&
     nrow(d_t_p>=nrow(d_tip)*0.5)){
    print(i)
    sel_trees[[as.character(n_sel_post[i])]]<-t1
    sel_trees_tips_lists[[as.character(n_sel_post[i])]]<-t1@phylo[["tip.label"]]
  }
}
sel_tree_tip<-unique(unlist(sel_trees_tips_lists))
name<-paste(args[3],"_",args[4],"_new",sep="")
Per_O_d<-data_clean%>%
  dplyr::filter(EPI_ISL%in%sel_tree_tip)
writeLines(sel_tree_tip,paste(name,".list",sep=""))
write.table(Per_O_d,paste(name,".tsv",sep=""),sep="\t",quote=F,row.names=F)
t2<-tree_subset(t_f_iq,node=as.numeric(names(sel_trees)[1]),levels_back=2)
Per_O_d2<-data_clean%>%
  dplyr::filter(EPI_ISL%in%t2@phylo[["tip.label"]])
writeLines(t2@phylo[["tip.label"]],paste(name,"2",".list",sep=""))
write.table(Per_O_d2,paste(name,"2",".tsv",sep=""),sep="\t",quote=F,row.names=F)
t1<-tree_subset(t_f_iq,node=as.numeric(names(sel_trees)[1]),levels_back=1)
Per_O_d1<-data_clean%>%
  dplyr::filter(EPI_ISL%in%t1@phylo[["tip.label"]])
writeLines(t1@phylo[["tip.label"]],paste(name,"3",".list",sep=""))
write.table(Per_O_d1,paste(name,"3",".tsv",sep=""),sep="\t",quote=F,row.names=F)
