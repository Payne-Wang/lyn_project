######################################################## 
#-------------------------------------------------------
# Topic: Microarray Data Analysis
# Author: Wang Haiquan
# Date: Thu Aug  1 17:14:10 2024
# Mail: mg1835020@smail.nju.edu.cn
#-------------------------------------------------------
########################################################

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("GEOquery", "limma", "pheatmap", "ggplot2", "EnhancedVolcano"))
# Load required packages
library(GEOquery)
library(limma)
library(pheatmap)
library(ggplot2)
library(EnhancedVolcano)
library(stringr)
library(dplyr)

# Download and load GSE167033 dataset
gse <- getGEO("GSE167033", GSEMatrix = TRUE, AnnotGPL = TRUE)
exprs_data <- exprs(gse[[1]])
pData <- pData(gse[[1]])
fData <- fData(gse[[1]])

# Quality control steps

# 1. Check for missing values
missing_values <- sum(is.na(exprs_data))
print(paste("Number of missing values:", missing_values))

# 2. Boxplot to view distribution across samples
boxplot(exprs_data, las = 2, main = "Distribution of expression values")

# 3. Density plot
plot(density(exprs_data[,1]), main = "Density plot of expression values", xlab = "Expression values")
for(i in 2:ncol(exprs_data)){
  lines(density(exprs_data[,i]), col = i)
}

# 4. Correlation heatmap
cor_matrix <- cor(exprs_data)
pheatmap(cor_matrix, main = "Sample Correlation Heatmap",
         annotation_col = pData["title"]%>%mutate(title=str_remove(title,"_rep[0-9]")))

# 5. PCA analysis
pca_result <- prcomp(t(exprs_data),scale. = T)
pca_data <- data.frame(PC1 = pca_result$x[,1], PC2 = pca_result$x[,2], group = pData$title%>%str_remove("_rep[0-9]"))
ggplot(pca_data, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size=3) +
  theme_bw() +
  ggtitle("PCA plot")

# Filter low expression genes
keep <- rowSums(exprs_data > 1) >= 5
exprs_data <- exprs_data[keep,]

# Normalization
exprs_data <- normalizeBetweenArrays(exprs_data)
boxplot(exprs_data, las = 2, main = "Distribution of expression values")

# Create design matrix
group <- factor(pData$`treatment:ch1`, levels = c("control", "treatment"))
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# Differential expression analysis
fit <- lmFit(exprs_data, design)
contrast.matrix <- makeContrasts(treatment - control, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Get differential expression gene results
results <- topTable(fit2, coef = 1, number = Inf)

# Visualization: Heatmap
top_genes <- head(rownames(results), 50)
pheatmap(exprs_data[top_genes,], scale = "row", show_rownames = FALSE,
         annotation_col = data.frame(Group = group))

# Visualization: Volcano plot
EnhancedVolcano(results,
                lab = rownames(results),
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 3.0,
                labSize = 6.0)

# Save results
write.csv(results, file = "GSE167033_differential_expression_results.csv")