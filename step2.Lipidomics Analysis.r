######################################################## 
#-------------------------------------------------------
# Topic: Lipidomics Analysis
# Author: Wang Haiquan
# Date: Sun Jul  7 19:25:01 2024
# Mail: mg1835020@smail.nju.edu.cn
#-------------------------------------------------------
########################################################
#devtools::install_github("junjunlab/bulkPseudotime")
library(limma)
library(bulkPseudotime)
library(dplyr)
library(stringr)
library(plyr)
library(tidyverse)
setwd("/data/whq/myproject/others_work/lyn_20240705/")
#-------------------------------------------------------
# Chapter:
#-------------------------------------------------------
# Read data and preprocess, using log method for standardization
lyn_tg_mat<-read.csv("SFF-CCl4-data_tidy-TG无e.csv")%>%
  column_to_rownames("LipidIon")%>%
  setNames(str_replace(colnames(.),"X","S"))%>%{log2(.+1)}

# Determine groups
lyn_tg_group<-read.csv("sampleList_lipid.csv")%>%
  mutate(samples=paste0("S",samples))%>%
  column_to_rownames("samples")%>%
  setNames("group")%>%
  mutate(group=factor(group,levels=c("Veh","12hrs","24hrs","48hrs")))

identical(colnames(lyn_tg_mat),rownames(lyn_tg_group))

#-------------------------------------------------------
# Chapter: Quality Control
#-------------------------------------------------------

# Quality Control

lyn_tg_mat%>%t()%>%prcomp(.,scale.=T) ->lyn_tg_pca
lyn_tg_pca%>%.$x%>%as.data.frame()%>%{
  ggplot()+
    geom_point(data = .,aes(x=PC1,y=PC2,color=lyn_tg_group$group))+
    geom_text_repel(data = .,aes(x=PC1,y=PC2),label=rownames(lyn_tg_group))+
    theme_prism()+
    labs(x=paste("PC1(",summary(lyn_tg_pca)%>%.$importance%>%.[2,1]*100,"%)",sep = ""),
         y=paste("PC2(",summary(lyn_tg_pca)%>%.$importance%>%.[2,2]*100,"%)",sep = ""))
}

pheatmap(cor(lyn_tg_mat),
         annotation_col = lyn_tg_group)

#-------------------------------------------------------
# Chapter: Calculate the relationship between 0, 12, 48h
#-------------------------------------------------------
lyn_tg_design<-model.matrix(~0+lyn_tg_group$group)
colnames(lyn_tg_design)<-paste0("a",levels(lyn_tg_group$group))
rownames(lyn_tg_design)<-colnames(lyn_tg_mat)

lyn_tg_contrast_mat_12vs48<-makeContrasts(a12hrs-a48hrs,levels=lyn_tg_design)
lyn_tg_contrast_mat_12vs0<-makeContrasts(a12hrs-aVeh,levels=lyn_tg_design)

lyn_tg_detg_12vs48<-lmFit(lyn_tg_mat,design = lyn_tg_design)%>%
  contrasts.fit(.,contrasts = lyn_tg_contrast_mat_12vs48)%>%
  eBayes()%>%
  topTable(number = Inf)%>%
  as.data.frame()%>%
  dplyr::select(logFC,P.Value,adj.P.Val)%>%
  setNames(c("log2FC","pvalue","padj"))%>%
  mutate(sig=ifelse(abs(log2FC)>1&pvalue<0.05,
                    ifelse(log2FC>1,"UP","DOWN"),"NO.Diff"),
         tg=rownames(.))


lyn_tg_detg_12vs0<-lmFit(lyn_tg_mat,design = lyn_tg_design)%>%
  contrasts.fit(.,contrasts = lyn_tg_contrast_mat_12vs0)%>%
  eBayes()%>%
  topTable(number = Inf)%>%
  as.data.frame()%>%
  dplyr::select(logFC,P.Value,adj.P.Val)%>%
  setNames(c("log2FC","pvalue","padj"))%>%
  mutate(sig=ifelse(abs(log2FC)>1&pvalue<0.05,
                    ifelse(log2FC>1,"UP","DOWN"),"NO.Diff"),
         tg=rownames(.))

ggplot(lyn_tg_detg_12vs48,aes(x=log2FC,y=-log10(pvalue),color=sig))+
  geom_point()+
  geom_hline(yintercept = -log10(0.05),lty=2,color="black")+
  geom_vline(xintercept = c(-1,1),lty=2,color="black")+
  scale_color_manual(values = c("UP"="red","DOWN"="blue","NO.Diff"="grey"))+
  scale_x_continuous(limits = c(-6,6))+
  theme_bw()+
  labs(title="12hrs vs 48hrs")

ggplot(lyn_tg_detg_12vs0,aes(x=log2FC,y=-log10(pvalue),color=sig))+
  geom_point()+
  geom_hline(yintercept = -log10(0.05),lty=2,color="black")+
  geom_vline(xintercept = c(-1,1),lty=2,color="black")+
  scale_color_manual(values = c("UP"="red","DOWN"="blue","NO.Diff"="grey"))+
  scale_x_continuous(limits = c(-6,6))+
  theme_bw()+
  labs(title="12hrs vs Veh")


lyn_tg_selected<-intersect(lyn_tg_detg_12vs48%>%filter(sig=="UP")%>%pull(tg),
                           lyn_tg_detg_12vs0%>%filter(sig=="UP")%>%pull(tg))
ggvenn::ggvenn(list("12hrs vs Veh (Up)"=lyn_tg_detg_12vs0%>%filter(sig=="UP")%>%pull(tg),
                    "12hrs vs 48hrs (Up)"=lyn_tg_detg_12vs48%>%filter(sig=="UP")%>%pull(tg)))



pheatmap::pheatmap(lyn_tg_mat[lyn_tg_selected,],scale = "row",
         cluster_cols = F,show_colnames = F,
         annotation_col = lyn_tg_group)

#-------------------------------------------------------
# Chapter: Analyze TG properties
#-------------------------------------------------------
lyn_tg_info<-rownames(lyn_tg_mat)%>%.[str_detect(.,"_")]%>%{
    a=.
    b=str_remove(a,"TG\$")%>%str_remove("\$")
    str_split(b,"_",simplify = T)%>%
      apply(.,1,function(x){str_split(x,":",simplify = T)%>%as.numeric})%>%t()%>%{rownames(.)=a;.}%>%
      as.data.frame()%>%
      setNames(c("FC","SC","TC","FDB","SDB","TDB"))
  }

lyn_tg_anno<-apply(lyn_tg_info,1,function(x){
  carbon_numbers=as.numeric(x[1:3])
  ifelse(all(carbon_numbers >= 6 & carbon_numbers <= 12),"MCT",
         ifelse(all(carbon_numbers > 12),"LCT",
                ifelse(any(carbon_numbers >= 6 & carbon_numbers <= 12) && any(carbon_numbers > 12),"MLCT","Mix")))
})%>%data.frame(row.names = names(.),
                TG=names(.),
                CT_type=.)
lyn_tg_anno$DU<-apply(lyn_tg_info,1,function(x){sum(x[4:6])})
lyn_tg_anno$DU_type<-apply(lyn_tg_info,1,function(x){
  double_bonds=x[4:6]
  saturated <- sum(double_bonds == 0)
  monounsaturated <- sum(double_bonds == 1)
  polyunsaturated <- sum(double_bonds > 1)
  if (saturated == 3) {
    saturation_type <- "STG"
  } else if (monounsaturated > 0 && polyunsaturated == 0) {
    saturation_type <- "MTG"
  } else if (polyunsaturated > 0) {
    saturation_type <- "PTG"
  } else {
    saturation_type <- "unknown_TG"
  }
  
  })

# Statistics of all TG molecules, pie chart
lyn_tg_anno%>%
  mutate(type=paste0(CT_type,"_",DU_type))%>%
  mutate(x="x")%>%{
    a=.
    #browser()
    ratio=ddply(a["type"],.(type),function(x){x%>%mutate(ratio=nrow(.)/nrow(a)*100,count=nrow(.))})%>%
      dplyr::select("type","ratio","count")%>%
      distinct()%>%
      arrange(ratio)
    b=merge(a,ratio,by="type",all.x=T)
    b$type=factor(b$type,levels=ratio$type,labels = paste0(ratio$type,"(",round(ratio$ratio,2),"%, count=",ratio$count,")"))
    b
  }%>%
  ggplot(aes(x=x,fill=type))+
  geom_bar(color="black")+
  coord_polar(theta = "y")+
  theme_bw()+
  theme(axis.title = element_blank(),axis.line = element_blank(),axis.text = element_blank(),panel.grid = element_blank(),axis.ticks = element_blank())

# Statistics of differential TG molecules, pie chart
lyn_tg_anno%>%
  mutate(type=paste0(CT_type,"_",DU_type))%>%
  filter(TG%in%lyn_tg_selected)%>%
  mutate(x="x")%>%{
    a=.
    #browser()
    ratio=ddply(a["type"],.(type),function(x){x%>%mutate(ratio=nrow(.)/nrow(a)*100,count=nrow(.))})%>%
      dplyr::select("type","ratio","count")%>%
      distinct()%>%
      arrange(ratio)
    b=merge(a,ratio,by="type",all.x=T)
    b$type=factor(b$type,levels=ratio$type,labels = paste0(ratio$type,"(",round(ratio$ratio,2),"%, count=",ratio$count,")"))
    b
  }%>%
  ggplot(aes(x=x,fill=type))+
  geom_bar(color="black")+
  coord_polar(theta = "y")+
  theme_bw()+
  theme(axis.title = element_blank(),axis.line = element_blank(),axis.text = element_blank(),panel.grid = element_blank(),axis.ticks = element_blank())


#-------------------------------------------------------
# Chapter: Calculate FFA chain content
#-------------------------------------------------------
lyn_tg_ffa<-lyn_tg_mat%>%.[rownames(lyn_tg_info),]%>%
  mutate(tg=rownames(.))%>%
  apply(.,1,function(x){
    list(x,x,x)%>%Reduce(rbind,.)%>%as.data.frame()%>%
      mutate(tg=str_remove(tg,"TG\$")%>%str_remove("\$")%>%str_split("_",simplify = F)%>%.[[1]])
  },simplify = F)%>%Reduce(rbind,.)%>%
  ddply(.,.(tg),function(x){
    x[,colnames(x)!="tg"]%>%as.matrix()%>%apply(.,2,function(y){sum(as.numeric(y))})
    
  })%>%column_to_rownames("tg")

pheatmap::pheatmap(lyn_tg_ffa,scale = "row",
                   cluster_cols = F,show_colnames = F,
                   annotation_col = lyn_tg_group)

# Screen through pseudotime
psetime_res <- my_bulkPseudotime(expMat = lyn_tg_ffa,group = lyn_tg_group$group,timePointCol = structure(rainbow(4),names=unique(lyn_tg_group$group)))
psetime_res@sample_pca_plot
set.seed(202412333)
psetime_res_heatmap<-pseudotime_heatmap(object = psetime_res,order = 3, heatmap_params = list(km = 3,
                                                                                              column_title="",
                                                                                              cluster_rows = T,
                                                                                              show_row_names = TRUE,
                                                                                              top_annotation = HeatmapAnnotation(timePoint = psetime_res@pseudotime_anno$group, 
                                                                                                                                 timeline = psetime_res@pseudotime_anno$time, border = T)))
psetime_res_heatmap = draw(psetime_res_heatmap)
psetime_res_tg_clust=row_dend(psetime_res_heatmap)
names(psetime_res_tg_clust)=paste0("C",names(psetime_res_tg_clust))
lapply(names(psetime_res_tg_clust),function(x){
  psetime_res_tg_clust[[x]]%>%
  as.hclust(.)->a
  data.frame(module=a$labels,cluster=x)
    
})%>%Reduce(rbind,.)->psetime_res_cluster


pseudotime_line(object = psetime_res,genes = psetime_res_cluster%>%filter(cluster=="C1")%>%pull(module))


# Relationship between saturation and length in each cluster

psetime_res_cluster%>%
  mutate(C_num=str_split(module,":",simplify = T)%>%.[,1]%>%as.numeric(),
         DU_num=str_split(module,":",simplify = T)%>%.[,2]%>%as.numeric())%>%
  ggplot(aes(x=C_num,y=DU_num,color=cluster,shape=cluster))+
  geom_point(size=5)+
  theme_bw()

#-------------------------------------------------------
# Chapter: Final processing and saving of results
#-------------------------------------------------------

# View variable sizes
require(dplyr)
require(stringr)
my_obj_in_env<-sapply(ls(),function(x){
  a=eval(expr = parse(text = x))%>%object.size()
  b=a/1024^3
  b%>%round(.,digits = 3)%>%as.numeric()})%>%
  data.frame(object=names(.),size_GB=.,row.names = NULL)%>%arrange(-size_GB)

my_obj_in_env[1:10,]

# Delete unnecessary variables
#rm()

# Save results
res_dir<-"/data/whq/myproject/others_work/lyn_20240705/RDS/"
image.name<-""
save.image(paste(res_dir,image.name,sep = ""))

lapply(ls(),function(x){filename=paste(res_dir,x,".rds",sep = "");saveRDS(eval(parse(text = x)),filename)})

# Custom function
my_bulkPseudotime<-function (expMat = NULL, timePointCol = NULL, timelineCol = NULL,group=NULL) 
{
  pca <- FactoMineR::PCA(t(expMat), scale.unit = T, ncp = 5, 
                         graph = F)
  res <- pca$ind$coord
  res <- res[, 1:2]
  dis <- stats::dist(res)
  dis <- as.matrix(dis)
  sample_dist_plot <- Heatmap(matrix = dis, cluster_rows = F, 
                              cluster_columns = F)
  res_p <- data.frame(res)
  res_p$sample <- rownames(res_p)
  raw.timeline <- c(0, lapply(1:ncol(expMat), function(x) {
    point <- 0 + dis[x - 1, x]
  }) %>% unlist() %>% cumsum())
  new.timeline <- scales::rescale(raw.timeline, to = c(0, 10))
  tpm.z <- as.data.frame(t(apply(expMat, 1, scale)))
  names(tpm.z) <- names(expMat)
  grid <- data.frame(time = seq(0, 10, length.out = 100 * ncol(expMat)))
  col_name <- colnames(expMat)
  cut_time <- cut(x = grid$time, breaks = new.timeline, labels = col_name[2:length(col_name)]) %>% 
    as.character()
  cut_time[1] <- col_name[1]
  grid_anno <- data.frame(time = grid$time, time_anno = cut_time)
  grid_anno$time_anno <- factor(grid_anno$time_anno, levels = colnames(expMat))
  grid_anno=merge(grid_anno,data.frame(time_anno=colnames(expMat),group=group),by="time_anno",all.x=T)%>%arrange(time)
  pseudotime.model.fun <- function(value) {
    time <- new.timeline
    data <- tibble::tibble(value = value, time = time)
    model <- stats::loess(value ~ time, data)
    predict <- grid %>% modelr::add_predictions(model)
    return(predict)
  }
  get_predicted_exp <- function(mat = NULL, timeline = NULL) {
    res <- apply(mat, 1, pseudotime.model.fun)
    results <- lapply(seq_along(res), function(x) {
      tmp <- res[[x]]$pred
      return(tmp)
    }) %>% do.call("rbind", .) %>% data.frame()
    rownames(results) <- row.names(mat)
    colnames(results) <- timeline$time
    return(results)
  }
  results <- get_predicted_exp(mat = tpm.z, timeline = grid)
  results_raw <- get_predicted_exp(mat = expMat, timeline = grid)
  if (is.null(timePointCol)) {
    time_col <- circlize::rand_color(n = length(unique(grid_anno$time_anno)))
  }
  else {
    time_col <- timePointCol
  }
  names(time_col) <- unique(grid_anno$group)
  if (is.null(timelineCol)) {
    col_fun = circlize::colorRamp2(c(0, 10), c("white", "#993399"))
  }
  else {
    col_fun = circlize::colorRamp2(seq(0, 10, length = length(timelineCol)), 
                                   timelineCol)
  }
  time_anno <- HeatmapAnnotation(timePoint = grid_anno$group, 
                                 timeline = grid_anno$time, col = list(timePoint = time_col, 
                                                                       timeline = col_fun), border = T)
  gene.pca <- FactoMineR::PCA(results, scale.unit = T, ncp = 5, 
                              graph = F)
  res2 <- gene.pca$ind$coord
  res2 <- res2[, 1:2]
  results2 <- cbind(results, res2)
  gene.pca2 <- FactoMineR::PCA(t(results_raw),  
                               ncp = 5, graph = F)
  res3 <- gene.pca2$ind$coord
  res3 <- data.frame(res3[, 1:2])
  res3$time_anno <- grid_anno$time_anno
  res3$group<-grid_anno$group
  results2_raw <- cbind(results_raw, res2) %>% tibble::rownames_to_column(var = "gene")
  results2_raw_long <- reshape2::melt(results2_raw, variable.name = "pseudotime", 
                                      value.name = "exp", id.vars = c("Dim.1", "Dim.2", "gene")) %>% 
    dplyr::mutate(pseudotime = as.numeric(as.character(pseudotime))) %>% 
    dplyr::left_join(y = grid_anno, by = c(pseudotime = "time"))
  
  sample_pca_plot <- ggplot() + 
    geom_text_repel(data = res_p, aes(x = Dim.1, y = Dim.2,label=sample))+
    geom_path(data = res_p, aes(x = Dim.1, y = Dim.2), linewidth = 1) + 
    geom_point(data = res_p, aes(x = Dim.1, y = Dim.2, color = group), size = 5) + 
    theme_bw(base_size = 12) + theme(panel.grid = element_blank()) 
  
  results2$atan2.1 <- atan2(results2$Dim.1, results2$Dim.2)
  results2$atan2.2 <- atan2(results2$Dim.1, -results2$Dim.2)
  results2$atan2.3 <- atan2(-results2$Dim.1, results2$Dim.2)
  results2$atan2.4 <- atan2(-results2$Dim.1, -results2$Dim.2)
  order1 <- dplyr::arrange(results2, results2$atan2.1)
  order2 <- dplyr::arrange(results2, results2$atan2.2)
  order3 <- dplyr::arrange(results2, results2$atan2.3)
  order4 <- dplyr::arrange(results2, results2$atan2.4)
  order_list <- list(order1, order2, order3, order4)
  names(order_list) <- paste0("Order", 1:4)
  ht_list <- lapply(1:4, function(x) {
    ht <- Heatmap(matrix = order_list[[x]][, 1:(100 * ncol(expMat))], 
                  name = "zscore", border = T, top_annotation = time_anno, 
                  cluster_rows = F, cluster_columns = F, show_row_names = F, 
                  show_column_names = F, column_title = paste0("Order", 
                                                               x))
    return(ht)
  })
  names(ht_list) <- paste0("Order", 1:4)
  ht_listcb <- Reduce("+", ht_list)
  res <- methods::new("bulkPseudotimeClass", sample_dist_plot = list(sample_dist_plot), 
                      sample_pca_plot = list(sample_pca_plot), pseudotime_matrix = results2, 
                      pseudotime_matrix_long = results2_raw_long, pseudotime_anno = grid_anno, 
                      ordered_matrix_list = list(order_list), heatmap_list = list(ht_list), 
                      heatmap_list_group = list(ht_listcb))
  return(res)
}