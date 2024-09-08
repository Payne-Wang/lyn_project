library(NMF)
library(ggplot2)
library(pheatmap)

# Tutorial: https://developer.aliyun.com/article/1256479
# Perform NMF decomposition
tg_data <- lyn_tg_mat
nmf_result_total <- nmf(tg_data, rank = 2:10, method = "brunet", nrun = 10, seed = 2024)
plot(nmf_result_total)
nmf_result <- nmf(tg_data, rank = 7, method = "brunet", nrun = 10, seed = 2024)

## Extract key genes
nmf_index <- extractFeatures(nmf_result, "max")

# Extract subtypes
nmf_predict <- predict(nmf_result)

# Extract results
W <- basis(nmf_result)
H <- coef(nmf_result)

# Convert decomposition results to data frames
W_df <- as.data.frame(W)
H_df <- as.data.frame(H)

# Add sample and TG names
rownames(W_df) <- rownames(tg_data)
colnames(H_df) <- colnames(tg_data)

# Plot heatmap to display NMF decomposition results
pheatmap::pheatmap(W_df, scale = "row", cluster_rows = TRUE, cluster_cols = TRUE,
                   main = "NMF Clustering of Samples",
                   color = colorRampPalette(c("navy", "white", "firebrick3"))(50))

pheatmap::pheatmap(H_df, scale = "column", cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = T, show_colnames = T,
                   main = "NMF Clustering of TG",
                   color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                   annotation_col = lyn_tg_group %>% mutate(nmf_predict = nmf_predict))

hdf_predict <- apply(H_df, 2, function(x) which.max(x))
wdf_predict <- apply(W_df, 1, function(x) which.max(x))

pheatmap(lyn_tg_mat,
         scale = "row",
         cluster_rows = T, cluster_cols = T,
         annotation_row = data.frame(row.names = rownames(lyn_tg_mat), group = paste0("C", wdf_predict)),
         annotation_col = lyn_tg_group %>% mutate(nmf_predict = paste0("C", nmf_predict)))

# Construct a pseudo-object
lyn_tg_nmf <- list(wide.res = lyn_tg_mat %>% mutate(gene = rownames(.), cluster = wdf_predict) %>% mutate(membership = 1),
                   long.res = lyn_tg_mat %>% mutate(gene = rownames(.), cluster = wdf_predict) %>% melt(id.vars = c("cluster", "gene")) %>%
                     setNames(c("cluster", "gene", "sample", "norm_value")) %>%
                     mutate(cluster_name = paste0("C", cluster)) %>%
                     merge(., lyn_tg_group %>% mutate(sample = rownames(.)), by = "sample", all.x = T) %>%
                     mutate(cell_type = group) %>% mutate(membership = 1),
                   type = "kmeans",
                   geneMode = "none",
                   geneType = "none")
visCluster(lyn_tg_nmf,
           plot.type = "line", ncol = 4)