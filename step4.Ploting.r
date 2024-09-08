######################################################## 
#-------------------------------------------------------
# Topic: Plotting
# Author: Wang Haiquan
# Date: Fri Jul 12 12:25:28 2024
# Mail: mg1835020@smail.nju.edu.cn
#-------------------------------------------------------
########################################################


# CCl4-0h: #f4aeb0
# CCl4-12h: #ff8080
# CCl4-24h: #a6cfce
# CCl4-48h: #65bab7

c("Veh"="#b3e2cd","12hrs"="#fdcdac","24hrs"="#cbd5e8","48hrs"="#f4cae4")
c("Veh"="#8EEAB5","12hrs"="#FFB77D","24hrs"="#A7B8F0","48hrs"="#FF9DE0")
c("Veh"="#66c2a5", "12hrs"="#fc8d62", "24hrs"="#8da0cb", "48hrs"="#e78ac3")


lyn_tg_detg_12vs48 <- readRDS("/data/whq/myproject/others_work/lyn_20240705/RDS/lyn_tg_detg_12vs48.rds")
lyn_tg_detg_12vs0 <- readRDS("/data/whq/myproject/others_work/lyn_20240705/RDS/lyn_tg_detg_12vs0.rds")
ggvenn::ggvenn(list("12hrs vs Veh (Up)"=lyn_tg_detg_12vs0%>%filter(sig=="UP")%>%pull(tg),
                    "12hrs vs 48hrs (Up)"=lyn_tg_detg_12vs48%>%filter(sig=="UP")%>%pull(tg)),fill_color = ggsci::pal_aaas()(2))

# Statistics of differential TG molecules, pie chart
lyn_tg_anno <- readRDS("/data/whq/myproject/others_work/lyn_20240705/RDS/lyn_tg_anno.rds")
lyn_tg_selected <- readRDS("/data/whq/myproject/others_work/lyn_20240705/RDS/lyn_tg_selected.rds")
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
  scale_fill_jama()+
  theme_bw()+
  theme(axis.title = element_blank(),axis.line = element_blank(),axis.text = element_blank(),panel.grid = element_blank(),axis.ticks = element_blank())




lyn_tg_pca <- readRDS("/data/whq/myproject/others_work/lyn_20240705/RDS/lyn_tg_pca.rds")
lyn_tg_group <- readRDS("/data/whq/myproject/others_work/lyn_20240705/RDS/lyn_tg_group.rds")


lyn_tg_pca%>%.$x%>%as.data.frame()%>%{
  ggplot()+
    geom_point(data = .,aes(x=PC1,y=PC2,fill=lyn_tg_group$group),size=5,pch=21,color="black")+
    geom_text_repel(data = .,aes(x=PC1,y=PC2),label=rownames(lyn_tg_group),max.overlaps = 5)+
    labs(fill="group")+
    scale_fill_manual(values =c("Veh"="#66c2a5", "12hrs"="#fc8d62", "24hrs"="#8da0cb", "48hrs"="#e78ac3"))+
    theme_bw()+
    labs(x=paste("PC1(",summary(lyn_tg_pca)%>%.$importance%>%.[2,1]*100,"%)",sep = ""),
         y=paste("PC2(",summary(lyn_tg_pca)%>%.$importance%>%.[2,2]*100,"%)",sep = ""))
    
}

lyn_tg_mat <- readRDS("/data/whq/myproject/others_work/lyn_20240705/RDS/lyn_tg_mat.rds")
lyn_tg_selected <- readRDS("/data/whq/myproject/others_work/lyn_20240705/RDS/lyn_tg_selected.rds")
pheatmap::pheatmap(lyn_tg_mat[lyn_tg_selected,],scale = "row",
                   cluster_cols = F,show_colnames = T,
                   color = colorRampPalette(c("skyblue","white","red"))(100),angle_col = 45,
                   annotation_col = lyn_tg_group,
                   annotation_colors = list(group=c("Veh"="#66c2a5", "12hrs"="#fc8d62", "24hrs"="#8da0cb", "48hrs"="#e78ac3")))



psetime_res <- readRDS("/data/whq/myproject/others_work/lyn_20240705/RDS/psetime_res.rds")
psetime_res@sample_pca_plot[[1]]+
  scale_color_manual(values =c("Veh"="#66c2a5", "12hrs"="#fc8d62", "24hrs"="#8da0cb", "48hrs"="#e78ac3"))

set.seed(202412333)
psetime_res_heatmap<-pseudotime_heatmap(object = psetime_res,order = 3, 
                                        heatmap_params = list(km = 3,
                                                              column_title="",
                                                              cluster_rows = T,
                                                              col = colorRampPalette(c("skyblue","white","red"))(100),
                                                              show_row_names = TRUE,
                                                              top_annotation = HeatmapAnnotation(timePoint = psetime_res@pseudotime_anno$group, 
                                                                                                 timeline = psetime_res@pseudotime_anno$time, border = T,
                                                                                                 col =list(timePoint=c("Veh"="#66c2a5", "12hrs"="#fc8d62", "24hrs"="#8da0cb", "48hrs"="#e78ac3") ))))

psetime_res_heatmap



# Relationship between saturation and length in each cluster
psetime_res_cluster <- readRDS("/data/whq/myproject/others_work/lyn_20240705/RDS/psetime_res_cluster.rds")
psetime_res_cluster%>%
  mutate(C_num=str_split(module,":",simplify = T)%>%.[,1]%>%as.numeric(),
         DU_num=str_split(module,":",simplify = T)%>%.[,2]%>%as.numeric())%>%
  ggplot(aes(x=C_num,y=DU_num,color=cluster,shape=cluster))+
  geom_point(size=5)+
  scale_color_npg()+
  theme_bw()