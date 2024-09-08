######################################################## 
#-------------------------------------------------------
# Topic: Further analysis of microarray data
# Author: Wang Haiquan
# Date: Fri Aug  2 13:51:26 2024
# Mail: mg1835020@smail.nju.edu.cn
#-------------------------------------------------------
########################################################

setwd("/data/whq/myproject/others_work/lyn_20240705/")
library(tidyverse)
library(ClusterGVis)
library(org.Mm.eg.db)
library(clusterProfiler)
library(openxlsx)
library(reshape2)
library(ggpubr)
library(Seurat)
library(ggrepel)
library(ggprism)

#-------------------------------------------------------
# Chapter: Load data
#-------------------------------------------------------
load("RDS/GSE167033.rdata")
identical(rownames(exprs_data),fData$ID)
exprs_data<-as.data.frame(exprs_data)
exprs_data$gene<-fData$`Gene symbol`
exprs_data<-exprs_data%>%filter(!duplicated(gene))%>%filter(!str_detect(gene,"\\/"))%>%remove_rownames%>%column_to_rownames(.,var = "gene")
lyn_array<-exprs_data
lyn_group<-pData%>%mutate(group=str_remove(title,"CCL4_")%>%str_remove("_rep[0-9]"))%>%
  mutate(group=factor(group,levels=c("ctrl","8h","1d","2d"),labels=paste0("CCl4-",c("0h","8h","24h","48h"))))

# Filter data, select only ctrl 8h 1d 2d
identical(colnames(lyn_array),rownames(lyn_group))
lyn_array<-lyn_array[,lyn_group$group%in%paste0("CCl4-",c("0h","8h","24h","48h"))]
lyn_group<-lyn_group[lyn_group$group%in%paste0("CCl4-",c("0h","8h","24h","48h")),]
identical(colnames(lyn_array),rownames(lyn_group))

lyn_array<-lyn_array[rowSums(lyn_array > 1) >= 10,]

# PCA
lyn_array%>%apply(.,1,scale)%>%prcomp(.) ->lyn_array_pca
lyn_array_pca%>%.$x%>%as.data.frame()%>%
  mutate(group=lyn_group$group,sample=paste0(lyn_group$group,"-",1:5))%>%
  column_to_rownames("sample")%>%{
    ggplot(.,aes(x=PC1,y=PC2,color=group))+
      geom_point(size=3)+
      geom_text_repel(aes(label=rownames(.)))+
      scale_color_manual(values = c("CCl4-0h"="#66c2a5", "CCl4-8h"="#fc8d62", "CCl4-24h"="#8da0cb", "CCl4-48h"="#e78ac3"))+
      theme_prism()+
      labs(x=paste("PC1(",summary(lyn_array_pca)%>%.$importance%>%.[2,1]*100,"%)",sep = ""),
           y=paste("PC2(",summary(lyn_array_pca)%>%.$importance%>%.[2,2]*100,"%)",sep = ""))+
      NoLegend()
  }

#-------------------------------------------------------
# Chapter: Further gene filtering
#-------------------------------------------------------
lyn_array_mean<-apply(lyn_array,1,function(x){
  tapply(x,lyn_group$group,mean)
})%>%t()%>%.[rowSums(.)>1,]%>%as.data.frame()%>%
  dplyr::select(paste0("CCl4-",c("0h","8h","24h","48h")))

lyn_array_mean_cm<-clusterData(exp = as.matrix(lyn_array_mean)%>%apply(.,1,scale)%>%t()%>%as.data.frame()%>%setNames(colnames(lyn_array_mean)),
                              cluster.method = "mfuzz",
                              cluster.num = 8,
                              seed = 2024,
                              scaleData = F)

# Remove genes with membership lower than 0.4
lyn_array_mean_cm$wide.res<-filter(lyn_array_mean_cm$wide.res,membership>0.4)
lyn_array_mean_cm$long.res<-filter(lyn_array_mean_cm$long.res,membership>0.4)%>%
   mutate(cluster_name=str_remove(cluster_name," \\(.{1,}\\)"))

visCluster(object = lyn_array_mean_cm,
           plot.type = "line",ncol = 4)

# Batch KEGG enrichment
lyn_array_mean_cm_enrichkegg<-lyn_array_mean_cm$wide.res%>%
  mutate(cluster=paste0("C",cluster))%>%
  mutate(cluster=factor(cluster,levels=paste0("C",1:8)))%>%
  split(.,.$cluster)%>%lapply(.,function(x){
    mapIds(org.Mm.eg.db,x$gene,"ENTREZID","SYMBOL")%>%
      enrichKEGG(.,"mmu","ncbi-geneid",pvalueCutoff=1)%>%
      setReadable(org.Mm.eg.db,"ENTREZID")
  })

# Select top 20 KEGG pathways
lyn_array_mean_cm_enrichkegg_df<-lyn_array_mean_cm_enrichkegg%>%
  lapply(.,function(x){
    x@result%>%.[1:30,]
  })%>%Reduce(rbind,.)%>%
  mutate(cluster=rep(names(lyn_array_mean_cm_enrichkegg),each=30))%>%
  dplyr::select("cluster",colnames(.)[colnames(.)!="cluster"])

# Save, download and highlight key information, display 5 terms for each cluster
# write.xlsx(lyn_array_mean_cm_enrichkegg_df,"lyn_array_mean_cm_enrichkegg_df_top10.xlsx")

# Construct candidate terms
lyn_array_mean_cm_enrichkegg_selected<-read.xlsx("lyn_array_mean_cm_enrichkegg_df_top30.xlsx",sheet = 2)%>%
  dplyr::select(cluster,Description)%>%
  mutate(Description=str_remove(Description," - Mus musculus \\(house mouse\\)"))%>%
  setNames(c("id","term"))

# Create heatmap
pdf('lyn_array_mean_cm_heatmap.pdf',height = 8,width = 10)
visCluster(object = lyn_array_mean_cm,
           plot.type = "both",
           show_row_dend = F,line.side = "left",
           column_names_rot = 0,
           annoTerm.data =lyn_array_mean_cm_enrichkegg_selected,
           ctAnno.col = ggsci::pal_aaas()(8),
           go.col = rep(ggsci::pal_aaas()(8),each=5),
           sample.col = c("CCl4-0h"="#66c2a5", "CCl4-8h"="#fc8d62", "CCl4-24h"="#8da0cb", "CCl4-48h"="#e78ac3")
)
dev.off()

#-------------------------------------------------------
# Chapter: Expression of some genes
#-------------------------------------------------------
lyn_array[c("Atf4","Eif2ak3","Sav1","Pik3ca","Braf","Mapk8"),]%>%
  mutate(gene=rownames(.))%>%
  melt(id.vars="gene")%>%
  merge(.,lyn_group%>%mutate(sample=rownames(.))%>%.[c("sample","group")],by.x="variable",by.y="sample",all.x=T)%>%
  setNames(c("sample","gene","expr","group"))%>%
  ggplot(aes(x=group,y=expr,fill=group))+
  geom_boxplot()+
  stat_compare_means(comparisons = list(c("CCl4-0h","CCl4-8h"),c("CCl4-8h","CCl4-24h"),c("CCl4-8h","CCl4-48h")),method = "t.test")+
  scale_fill_manual(values = c("CCl4-0h"="#66c2a5", "CCl4-8h"="#fc8d62", "CCl4-24h"="#8da0cb", "CCl4-48h"="#e78ac3"))+
  facet_wrap(~gene,scales = "free",ncol = 3)+
  theme_bw()+NoLegend()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))+
  labs(x="",y="log2(expr+1)")

lyn_array[c("Ccna2","Cdk2","Pcna","Mcm2"),]%>%
  mutate(gene=rownames(.))%>%
  melt(id.vars="gene")%>%
  merge(.,lyn_group%>%mutate(sample=rownames(.))%>%.[c("sample","group")],by.x="variable",by.y="sample",all.x=T)%>%
  setNames(c("sample","gene","expr","group"))%>%
  ggplot(aes(x=group,y=expr,fill=group))+
  geom_boxplot()+
  stat_compare_means(comparisons = list(c("CCl4-0h","CCl4-8h"),c("CCl4-8h","CCl4-24h"),c("CCl4-8h","CCl4-48h")),method = "t.test")+
  scale_fill_manual(values = c("CCl4-0h"="#66c2a5", "CCl4-8h"="#fc8d62", "CCl4-24h"="#8da0cb", "CCl4-48h"="#e78ac3"))+
  facet_wrap(~gene,scales = "free",ncol = 2)+
  theme_bw()+NoLegend()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))+
  labs(x="",y="log2(expr+1)")

lyn_array[c("Acaca","Acly","Acss2","Fasn"),]%>%
  mutate(gene=rownames(.))%>%
  melt(id.vars="gene")%>%
  merge(.,lyn_group%>%mutate(sample=rownames(.))%>%.[c("sample","group")],by.x="variable",by.y="sample",all.x=T)%>%
  setNames(c("sample","gene","expr","group"))%>%
  ggplot(aes(x=group,y=expr,fill=group))+
  geom_boxplot()+
  stat_compare_means(comparisons = list(c("CCl4-0h","CCl4-8h"),c("CCl4-8h","CCl4-24h"),c("CCl4-8h","CCl4-48h")),method = "t.test")+
  scale_fill_manual(values = c("CCl4-0h"="#66c2a5", "CCl4-8h"="#fc8d62", "CCl4-24h"="#8da0cb", "CCl4-48h"="#e78ac3"))+
  facet_wrap(~gene,scales = "free",ncol = 2)+
  theme_bw()+NoLegend()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))+
  labs(x="",y="log2(expr+1)")

#-------------------------------------------------------
# Chapter: Final processing and saving of results
#-------------------------------------------------------

# Check variable sizes
require(dplyr)
require(stringr)
my_obj_in_env<-sapply(ls(),function(x){
  a=eval(expr = parse(text = x))%>%object.size()
  b=a/1024^3
  b%>%round(.,digits = 3)%>%as.numeric()})%>%
  data.frame(object=names(.),size_GB=.,row.names = NULL)%>%arrange(-size_GB)

my_obj_in_env[1:10,]

# Remove unnecessary variables
# rm()

# Save results
res_dir<-"/data/whq/myproject/others_work/lyn_20240705/RDS/"
image.name<-""
save.image(paste(res_dir,image.name,sep = ""))

lapply(ls(),function(x){filename=paste(res_dir,x,".rds",sep = "");saveRDS(eval(parse(text = x)),filename)})