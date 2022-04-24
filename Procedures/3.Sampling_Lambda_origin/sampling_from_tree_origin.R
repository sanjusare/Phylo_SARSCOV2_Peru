args<-commandArgs(trailingOnly = TRUE)
Sys.setlocale("LC_TIME", "English")
if(length(args)>2| length(args)<2){
  print("usage=Rscript sampling_from_tree_origin.R [metadata_file] [tree_file]")
  break
}
if (!require(dplyr, quietly = T)) install.packages("dplyr")
library(dplyr)
if (!require(tidyr, quietly = T)) install.packages("tidyr")
library(tidyr)
data<-read.table(args[1],sep="\t",quote="",header=TRUE)
data$name_analysis<-gsub("'","_",data$name_analysis)
data$name_analysis<-gsub("Coted_Ivoire","CotedIvoire",data$name_analysis)
data$C37<-"no"
data$C37[grepl("Lambda",data$variant)]<-"yes"
data$date<-as.Date(data$date)
t<-readLines(args[2])
t_cl<-gsub("'","",t)
writeLines(t_cl,"test.treefile")
t_f_iq<-read.iqtree("test.treefile")
data_clean<-data%>%
  dplyr::filter(name_analysis%in%t_f_iq@phylo[["tip.label"]])
p_t_f_iq<-ggtree(t_f_iq)%<+%data_clean
nodes<-p_t_f_iq[["data"]]%>%
  dplyr::filter(C37=="yes")
r_node<-p_t_f_iq[["data"]]$parent[p_t_f_iq[["data"]]$node==p_t_f_iq[["data"]]$parent]
#create the pathway (parents) of each of the Lambda tips
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
lam_n_pat<-Reduce(intersect,saved_nodes_list)
lam_n<-max(lam_n_pat)

n_sel_post<-c()
for(i in lam_n_pat){
  if(i==r_node){
  }else if(p_t_f_iq[["data"]]$SH_aLRT[p_t_f_iq[["data"]]$node==i]>=0.9){
    n_sel_post<-append(n_sel_post,i)
  }
}

for(i in 1:length(n_sel_post)){
  t1<-tree_subset(t_f_iq,node=n_sel_post[i],levels_back=0)
  n_tip<-length(t1@phylo[["tip.label"]])
  if(n_tip>=200){
    sel_tree<-t1
    sel_trees_tip<-t1@phylo[["tip.label"]]
    break
  }
}

rest<-200-nrow(nodes)

#create the metadata with the tips
filt_data<-data_clean%>%
  dplyr::filter(name_analysis%in%sel_trees_tip)
#create a metadata with lambda genomes
C37_fil_d<-filt_data%>%
  filter(C37=="yes")
#create a metadata with all the peruvian genomes different to lambda
Peru_fil_d<-filt_data%>%
  filter(C37=="no",Country=="Peru")
#create the data with the number of genomes to extract from the no lambda peruvian data
nc37_Peru_g<-Peru_fil_d%>%
  group_by(pango_lineage)%>%
  summarise(count=n())%>%
  mutate(percent=count/sum(count))%>%
  mutate(sample_n=round(percent*(rest/2),0))
nc37_Peru_g$sample_n[nc37_Peru_g$sample_n==0]<-1
#sampling of the genomes no lambda peru
sample_nc37_Peru_list<-list()
for (i in 1:nrow(nc37_Peru_g)){
  row<-nc37_Peru_g[i,] #crea obj de 1 fila
  sample_nc37_Peru_list[[row$pango_lineage]]<-Peru_fil_d%>%
    dplyr::filter(pango_lineage==row$pango_lineage)%>%
    dplyr::sample_n(size=row$sample_n)
}
sample_nc37_Peru<-do.call(rbind,sample_nc37_Peru_list)
#create a metadata with no lambda, no peru genomes
nC37_fil_d<-filt_data%>%
  filter(C37=="no",Country!="Peru")
#create the data with the number of genomes to extract from the no lambda no peru data
n_c37_g<-nC37_fil_d%>%
  group_by(pango_lineage)%>%
  summarise(count=n())%>%
  mutate(percent=count/sum(count))%>%
  mutate(sample_n=round(percent*(rest-nrow(sample_nc37_Peru)),0))
n_c37_g$sample_n[n_c37_g$sample_n==0]<-1
#sampling of the genomes no lambda no peru
sample_nC37_list<-list()
for (i in 1:nrow(n_c37_g)){
  row<-n_c37_g[i,] #crea obj de 1 fila
  sample_nC37_list[[row$pango_lineage]]<-nC37_fil_d%>%
    dplyr::filter(pango_lineage==row$pango_lineage)%>%
    dplyr::sample_n(size=row$sample_n)
}
sample_nC37<-do.call(rbind,sample_nC37_list)
#join the samples
final_sample_prev<-rbind(C37_fil_d,sample_nC37)
final_sample<-rbind(final_sample_prev,sample_nc37_Peru)
#add the information of location
final_sample$location<-"World"
final_sample$location[final_sample$Country=="Peru"]<-"Peru"
#save the important documents
final_sample_sel<-final_sample%>%
  select(name_analysis,location)
name<-paste(sapply(strsplit(args[1],"\\."),"[",1),"_treesample",sep="")
writeLines(final_sample$name_analysis,paste(name,".list",sep=""))
write.table(final_sample_sel,paste(name,".tsv",sep=""),sep="\t",quote=F,row.names=F)
