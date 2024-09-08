######################################################## 
#-------------------------------------------------------
# Topic: Li Yining Data Analysis -- Transcriptome + Metabolome
# Author: Wang Haiquan
# Date: Fri Jul  5 14:36:14 2024
# Mail: mg1835020@smail.nju.edu.cn
#-------------------------------------------------------
########################################################
library(TCseq)
library(ClusterGVis)
library(Biobase)
library(ggsci)
library(openxlsx)
#devtools::install_github("junjunlab/ClusterGVis")
setwd("/data/whq/myproject/others_work/lyn_20240705/")
#-------------------------------------------------------
# Chapter: Data Input
#-------------------------------------------------------
# Read data

lyn_count <- read.xlsx("RNAseq_countsandFPKM.xlsx") %>% filter(!duplicated(GeneSymbol)) %>%
  column_to_rownames("GeneSymbol") %>%
  dplyr::select(colnames(.) %>% .[str_detect(.,"^(v|([0-9]{1,}(h|d)))-[1-3]$")])

lyn_rpkm <- read.xlsx("RNAseq_countsandFPKM.xlsx") %>% filter(!duplicated(GeneSymbol)) %>%
  column_to_rownames("GeneSymbol") %>%
  dplyr::select(colnames(.) %>% .[str_detect(.,"^(v|([0-9]{1,}(h|d)))-[1-3]_FPKM")]) %>%
  setNames(colnames(.) %>% str_remove("_FPKM"))

# Define groups
lyn_rna_group <- data.frame(row.names = colnames(lyn_count),
                          group=colnames(lyn_count) %>% str_remove("-[0-3]") %>% factor(.,levels=c("v","12h","1d","2d"),label=c("Veh","12hrs","24hrs","48hrs"))) %>%
  .[order(.$group),,drop=F]

# Rearrange data, log transform RPKM
lyn_count <- lyn_count[,rownames(lyn_rna_group)] %>% .[rowSums(lyn_count)>5,]
lyn_rpkm <- lyn_rpkm[,rownames(lyn_rna_group)] %>% .[rownames(lyn_count),] %>% {log2(.+1)}
identical(colnames(lyn_rpkm),colnames(lyn_count))
identical(rownames(lyn_rpkm),rownames(lyn_count))
#-------------------------------------------------------
# Chapter: Quality Control
#-------------------------------------------------------
# Create DGElist object

lyn_DGElist <- lyn_count %>% DGEList(.,group = lyn_rna_group$group) %>% calcNormFactors()

# Quality control

lyn_rpkm %>% t() %>% prcomp(.,scale.=T) -> lyn_rpkm_pca
lyn_rpkm_pca %>% .$x %>% as.data.frame() %>%
  ggplot(aes(x=PC1,y=PC2,color=lyn_rna_group$group))+
  geom_point()+
  geom_text_repel(aes(label=rownames(lyn_rna_group)))+
  theme_prism()+
  labs(x=paste("PC1(",summary(lyn_rpkm_pca) %>% .$importance %>% .[2,1]*100,"%)",sep = ""),
       y=paste("PC2(",summary(lyn_rpkm_pca) %>% .$importance %>% .[2,2]*100,"%)",sep = ""))

pheatmap(cor(lyn_rpkm),
         annotation_col = lyn_rna_group)

#-------------------------------------------------------
# Chapter: Trend Analysis
#-------------------------------------------------------

# Further filter genes
lyn_rpkm_mean <- apply(lyn_rpkm,1,function(x){
  tapply(x,lyn_rna_group$group,mean)
}) %>% t() %>% .[rowSums(.)>1,]

#getClusters(exp = lyn_rpkm_mean)

lyn_rpkm_mean_cm <- clusterData(exp = as.matrix(lyn_rpkm_mean),
                              cluster.method = "mfuzz",
                              cluster.num = 8,
                              seed = 2024,
                              scaleData = T)

# Remove genes with membership lower than 0.4
lyn_rpkm_mean_cm$wide.res <- filter(lyn_rpkm_mean_cm$wide.res,membership>0.4)
lyn_rpkm_mean_cm$long.res <- filter(lyn_rpkm_mean_cm$long.res,membership>0.4) %>%
  mutate(cluster_name=str_remove(cluster_name," \\(.{1,}\\)"))

visCluster(object = lyn_rpkm_mean_cm,
           plot.type = "line",ncol = 4)


# Batch KEGG enrichment
lyn_rpkm_mean_cm_enrichkegg <- lyn_rpkm_mean_cm$wide.res %>%
  mutate(cluster=paste0("C",cluster)) %>%
  mutate(cluster=factor(cluster,levels=paste0("C",1:8))) %>%
  split(.,.$cluster) %>% lapply(.,function(x){
  mapIds(org.Mm.eg.db,x$gene,"ENTREZID","SYMBOL") %>%
    enrichKEGG(.,"mmu","ncbi-geneid",pvalueCutoff=1) %>%
    setReadable(org.Mm.eg.db,"ENTREZID")
})


# Select top 20 KEGG pathways
lyn_rpkm_mean_cm_enrichkegg_df <- lyn_rpkm_mean_cm_enrichkegg %>%
  lapply(.,function(x){
    x@result %>% .[1:20,]
  }) %>% Reduce(rbind,.) %>%
  mutate(cluster=rep(names(lyn_rpkm_mean_cm_enrichkegg),each=20)) %>%
  dplyr::select("cluster",colnames(.)[colnames(.)!="cluster"])

# Save, download and highlight key information, select 5 terms to display for each cluster
#write.xlsx(lyn_rpkm_mean_cm_enrichkegg_df,"lyn_rpkm_mean_cm_enrichkegg_df_top10.xlsx")

# Construct candidate terms
lyn_rpkm_mean_cm_enrichkegg_selected <- read.xlsx("lyn_rpkm_mean_cm_enrichkegg_df_top10.xlsx",sheet = 3) %>%
  dplyr::select(cluster,Description) %>%
  mutate(Description=str_remove(Description," - Mus musculus \\(house mouse\\)")) %>%
  setNames(c("id","term"))

# Create heatmap
pdf('lyn_rpkm_mean_cm_heatmap.pdf',height = 8,width = 10)
visCluster(object = lyn_rpkm_mean_cm,
           plot.type = "both",
           show_row_dend = F,line.side = "left",
           column_names_rot = 0,
           annoTerm.data =lyn_rpkm_mean_cm_enrichkegg_selected,
           ctAnno.col = ggsci::pal_aaas()(8),
           go.col = rep(ggsci::pal_aaas()(8),each=5),
           sample.col = c("Veh"="#66c2a5", "12hrs"="#fc8d62", "24hrs"="#8da0cb", "48hrs"="#e78ac3")
           )
dev.off()

#-------------------------------------------------------
# Chapter: Expression of Some Genes
#-------------------------------------------------------
lyn_rpkm[c("Dvl3","Fzd9","Wnt9a","Yap1",
           "Gsk3b","Smad2","Smad3","Smad4"),] %>%
  mutate(gene=rownames(.)) %>%
  melt(id.vars="gene") %>%
  merge(.,lyn_rna_group %>% mutate(sample=rownames(.)),by.x="variable",by.y="sample",all.x=T) %>%
  setNames(c("sample","gene","rpkm","group")) %>%
  ggplot(aes(x=group,y=rpkm,fill=group))+
  geom_boxplot()+
  stat_compare_means(comparisons = list(c("Veh","12hrs"),c("12hrs","24hrs"),c("12hrs","48hrs")),method = "t.test")+
  scale_fill_manual(values = c("Veh"="#66c2a5", "12hrs"="#fc8d62", "24hrs"="#8da0cb", "48hrs"="#e78ac3"))+
  facet_wrap(~gene,scales = "free",ncol = 4)+
  theme_bw()+NoLegend()+
  labs(x="",y="log2(FPKM+1)")


lyn_rpkm[c("Acaca","Acly","Acss2","Fasn"),] %>%
  mutate(gene=rownames(.)) %>%
  melt(id.vars="gene") %>%
  merge(.,lyn_rna_group %>% mutate(sample=rownames(.)),by.x="variable",by.y="sample",all.x=T) %>%
  setNames(c("sample","gene","rpkm","group")) %>%
  ggplot(aes(x=group,y=rpkm,fill=group))+
  geom_boxplot()+
  stat_compare_means(comparisons = list(c("Veh","12hrs"),c("Veh","24hrs"),c("Veh","48hrs")),method = "t.test")+
  scale_fill_manual(values = c("Veh"="#66c2a5", "12hrs"="#fc8d62", "24hrs"="#8da0cb", "48hrs"="#e78ac3"))+
  facet_wrap(~gene,scales = "free",ncol = 4)+
  theme_bw()+NoLegend()+
  labs(x="",y="log2(FPKM+1)")

#-------------------------------------------------------
# Chapter: Final Processing and Saving of Results
#-------------------------------------------------------

# Check variable sizes
require(dplyr)
require(stringr)
my_obj_in_env <- sapply(ls(),function(x){
  a=eval(expr = parse(text = x)) %>% object.size()
  b=a/1024^3
  b %>% round(.,digits = 3) %>% as.numeric()}) %>%
  data.frame(object=names(.),size_GB=.,row.names = NULL) %>% arrange(-size_GB)

my_obj_in_env[1:10,]

# Delete unnecessary variables
#rm()

# Save results
res_dir <- "/data/whq/myproject/others_work/lyn_20240705/RDS/"
image.name <- ""
save.image(paste(res_dir,image.name,sep = ""))

lapply(ls(),function(x){filename=paste(res_dir,x,".rds",sep = "");saveRDS(eval(parse(text = x)),filename)})