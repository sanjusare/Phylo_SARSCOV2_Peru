#1:set working directory####
args<-commandArgs(trailingOnly = TRUE)
#2:load libraries and establish settings####
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
theme_set(theme_bw(base_size=12))
Sys.setlocale("LC_TIME", "English")
colreg<-c("center"="firebrick","north"="dodgerblue","south"="forestgreen","mid-east"="darkmagenta","north-east"="goldenrod","south-east"="hotpink3")
#3:read files####
#read the epidemiological data (from section 3_4, step 5)
f_df_cas_all<-read.csv(file=args[1],sep="\t",header=TRUE)%>%
  dplyr::mutate(collection_date=as.Date(collection_date))
#read the genome data with sublineages information (from section 3_4, step 5)
df_p_reg<-read.csv(file=args[2],sep="\t",header=TRUE)
#4:plot of correlation between Gamma sublineages genomes and estimated cases by week by region####
p_cor_leg<-ggscatter(subset(f_df_cas_all,count!=0),x="count",y="smooth_lin_cases",
                     color="black",fill="Region",shape=21,size=3,
                     add="reg.line", 
                     add.params=list(color="blue",fill="lightgray"),
                     conf.int=TRUE,cor.coef=TRUE,
                     cor.coeff.args=list(method="pearson",label.x=3,label.sep="\n"))+
  labs(y="Number of Cases",x="Number of genomes")+
  scale_fill_manual(values=colreg)+
  facet_wrap(~SubGs,scales="free")+
  theme(legend.position="bottom")+
  guides(fill=guide_legend(nrow=1))
cor_leg<-get_legend(p_cor_leg)
p_cor<-p_cor_leg+
  theme(legend.position="none")
#5:sampling by sublineage according to the number of cases by region by week####
sel_gens_list<-list()
gen_group_cases_list<-list()
for(i in levels(as.factor(f_df_cas_all$SubGs))){
  print(i)
  f_df_cas_all_t<-f_df_cas_all%>%
    dplyr::filter(SubGs==i)%>%
    dplyr::ungroup()%>%
    dplyr::mutate(samp_gen=round(smooth_lin_cases/(sum(smooth_lin_cases,na.rm=T)/200),0),
                  collection_date=as.Date(collection_date))
  f_df_cas_all_t$samp_gen[is.na(f_df_cas_all_t$samp_gen)]<-0
  datos_list<-list()
  for(z in 1:nrow(f_df_cas_all_t)){
    print(z)
    set.seed(234364151)
    if(f_df_cas_all_t$samp_gen[z]>0){
      name<-paste(i,f_df_cas_all_t$Region[z],f_df_cas_all_t$collection_date[z],sep="_")
      te<-df_p_reg%>%
        dplyr::mutate(collection_date=as.Date(collection_date))%>%
        dplyr::filter(Region==f_df_cas_all_t$Region[z]&
                        collection_date==f_df_cas_all_t$collection_date[z]&
                        SubGs==i)
      if(nrow(te)>=f_df_cas_all_t$samp_gen[z]){
        datos_list[[name]]<-te%>%
          dplyr::sample_n(f_df_cas_all_t$samp_gen[z])
      }else{
        datos_list[[name]]<-te
      }
    }
  }
  sel_gens_list[[i]]<-do.call(rbind,datos_list)
  gen_group_cases_list[[i]]<-sel_gens_list[[i]]%>%
    dplyr::group_by(Region,collection_date,SubGs)%>%
    dplyr::summarise(samp_count=n())%>%
    dplyr::right_join(f_df_cas_all_t)
  gen_group_cases_list[[i]]$samp_count[is.na(gen_group_cases_list[[i]]$samp_count)]<-0
}
gen_group_cases<-do.call(rbind,gen_group_cases_list)
#6:plot of correlation between sampled Gamma sublineages genomes and estimated cases by week by region after sampling####
p_cor_sam<-ggscatter(subset(gen_group_cases,count!=0),x="samp_count",y="smooth_lin_cases",
                     color="black",fill="Region",shape=21,size=3,
                     add="reg.line", 
                     add.params=list(color="blue",fill="lightgray"),
                     conf.int=TRUE,cor.coef=TRUE,
                     cor.coeff.args=list(method="pearson",label.x=1,label.sep="\n"))+
  labs(y="Number of Cases",x="Number of genomes")+
  scale_fill_manual(values=colreg)+
  facet_wrap(~SubGs,scales="free")+
  theme(legend.position="none")
p_cor_fin<-ggarrange(p_cor,p_cor_sam,labels=c("A","B"),nrow=2,ncol=1,common.legend=T,legend="bottom",legend.grob=cor_leg)
ggsave("Gamma_sampling_correlation.jpg",plot=p_cor_fin,dpi=500,height=10,width=7.5)
#7:Plot to evaluate cases, available genomes and sampled genomes####
for(i in levels(as.factor(gen_group_cases$SubGs))){
  d_plot<-gen_group_cases%>%
    dplyr::filter(SubGs==i)
  factor<-max(gen_group_cases$count,na.rm=T)
  p_sampled<-ggplot(d_plot,aes(x=as.POSIXct(collection_date)))+
    geom_bar(aes(y=count),stat="identity",fill="lightgrey")+
    geom_bar(aes(y=samp_count,fill=Region),stat="identity")+
    geom_bar(data=subset(d_plot,samp_gen>0),aes(y=samp_gen),stat="identity",color="red",fill="transparent",linewidth=0.2)+
    geom_line(aes(y=smooth_lin_cases/factor))+
    facet_wrap(~Region,nrow=3,ncol=2,scales="free_y")+
    scale_y_continuous(sec.axis=sec_axis(trans=~.*factor,name=paste(i,"cases")))+
    scale_x_datetime(date_breaks="1 month",date_labels="%Y-%b",limits=c(as.POSIXct("2020-10-01"),as.POSIXct("2022-01-31")))+
    theme(legend.position="none",
          axis.title.x=element_blank(),
          axis.text.x=element_text(angle=90,hjust=0,vjust=0.5))+
    scale_fill_manual(values=colreg)+
    labs(y="Number of genomes")
  ggsave(paste(i,"_sampling_distribution.jpg",sep=""),plot=p_sampled,dpi=300,height=7.5,width=7.5)
}
#8:Saving the sampled genomes and tables with relevant information for analysis####
for(i in 1:length(sel_gens_list)){
  print(names(sel_gens_list[i]))
  lin<-names(sel_gens_list[i])
  d_reg<-sel_gens_list[[i]]%>%
    dplyr::select(name_analysis,Region)
  writeLines(sel_gens_list[[i]]$name_analysis,paste(lin,"_BI_fin_v2.list",sep=""))
  write.table(sel_gens_list[[i]],paste(lin,"_meta_fin_v2.tsv",sep=""),sep="\t",quote=F,row.names=F)
  write.table(d_reg,paste(lin,"_region_fin_v2.tsv",sep=""),sep="\t",quote=F,row.names=F)
}