args<-commandArgs(trailingOnly = TRUE)
Sys.setlocale("LC_TIME", "English")
if(length(args)>2| length(args)<2){
  print("usage=Rscript sampling_from_tree_introduction.R [metadata_file] [tree_file]")
  break
}
if (!require(dplyr, quietly = T)) install.packages("dplyr")
library(dplyr)
if (!require(tidyr, quietly = T)) install.packages("tidyr")
library(tidyr)
if (!require(ggtree, quietly = T)) BiocManager::install("ggtree")
library(ggtree)
data<-read.table(args[1],sep="\t",quote="",header=TRUE)
data$name_analysis<-gsub("'","_",data$name_analysis)
data$Peru<-"no"
data$Peru[data$Country=="Peru"]<-"yes"
data$date<-as.Date(data$date)
t_f_iq<-read.iqtree(args[2])
data_clean<-data%>%
  dplyr::filter(name_analysis%in%t_f_iq@phylo[["tip.label"]])
p_t_f_iq<-ggtree(t_f_iq)%<+%data_clean
nodes<-p_t_f_iq[["data"]]%>%
  dplyr::filter(Country=="Peru")
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
n_sel_post<-c()
p_t_f_iq[["data"]]$SH_aLRT[is.na(p_t_f_iq[["data"]]$SH_aLRT)]<-0
for(i in saved_nodes_unique){
  if(i==r_node){
  }else if((p_t_f_iq[["data"]]$SH_aLRT[p_t_f_iq[["data"]]$node==i]*100)>=70){
    n_sel_post<-append(n_sel_post,i)
  }
}
sel_trees<-list()
sel_trees_tips_lists<-list()
sel_trees_big<-list()
for(i in 1:length(n_sel_post)){
  t1<-tree_subset(t_f_iq,node=n_sel_post[i],levels_back=0)
  n_tip<-length(t1@phylo[["tip.label"]])
  if(n_tip<=20){
    sel_trees[[as.character(n_sel_post[i])]]<-t1
    sel_trees_tips_lists[[as.character(n_sel_post[i])]]<-t1@phylo[["tip.label"]]
  }else{
    sel_trees_big[[as.character(n_sel_post[i])]]<-t1
  }
}
sel_tree_tip<-unique(unlist(sel_trees_tips_lists))
abs_nodes<-nodes%>%
  dplyr::filter(!label%in%sel_tree_tip)
rest<-200-length(sel_tree_tip)-nrow(abs_nodes)
if(nrow(abs_nodes>0)){
#anotar los nodos bien soportados de clados grandes y muestrear de ahi lo faltante  
  big_n_sam<-list()
  for(i in 1:nrow(abs_nodes)){
    new_nodes<-saved_nodes_list[[as.character(abs_nodes$node[i])]]
    new_sel_nodes<-new_nodes[new_nodes%in%n_sel_post]
    new_f_sel_nodes<-new_sel_nodes[!(new_sel_nodes%in%as.numeric(names(sel_trees)))]
    big_n_sam[as.character(abs_nodes$node[i])]<-new_f_sel_nodes[1]
  }
  big_n_sam_un<-unique(unlist(big_n_sam))
  sam_size<-round(rest/length(big_n_sam_un),0)
  big_tips_list<-list()
  for(i in 1:length(big_n_sam_un)){
    big_tips_list[[as.character(big_n_sam_un[i])]]<-p_t_f_iq[["data"]]%>%
      dplyr::filter(isTip==TRUE)%>%
      dplyr::filter(!(label%in%sel_tree_tip))%>%
      dplyr::filter(label%in%sel_trees_big[[as.character(big_n_sam_un[i])]]@phylo[["tip.label"]])%>%
      dplyr::filter(Country!="Peru")%>%
      dplyr::sample_n(sam_size)
  }
  big_tips<-do.call(rbind,big_tips_list)
  big_tips_f<-big_tips$label
}else{
#seleccionar al azar genomas faltantes que no hayan estado previamente muestreados
  big_tips<-p_t_f_iq[["data"]]%>%
    dplyr::filter(isTip==TRUE)%>%
    dplyr::filter(!(label%in%sel_tree_tip))%>%
    dplyr::filter(Country!="Peru")%>%
    dplyr::sample_n(rest)
  big_tips_f<-big_tips$label
}
final_tips<-c(big_tips_f,sel_tree_tip,abs_nodes$label)
data_clean_fn<-data_clean%>%
  dplyr::filter(name_analysis%in%final_tips)
name<-paste(sapply(strsplit(args[1],"\\."),"[",1),"_treesample",sep="")
writeLines(final_tips,paste(name,".list",sep=""))
write.table(data_clean_fn,paste(name,".tsv",sep=""),sep="\t",quote=F,row.names=F)
