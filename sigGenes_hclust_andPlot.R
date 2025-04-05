library(stats)
library(tibble)
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)
library(data.table)

data <- readRDS("/scratch/labs/khatrilab/aleewq/data/fresh_3m052_ciceroobj_logGA_allcells.RDS")
specific_cell_subset <- subset(data, idents = c("ly6c+, Day 28"))
specific_countMat <- GetAssayData(specific_cell_subset, slot = "counts")
cell_subset_d1 <- subset(data, idents = c("ly6c+, Day 1"))
countMat_d1 <- GetAssayData(cell_subset_d1, slot = "counts")
cell_subset_d0 <- subset(data, idents = c("ly6c+, Day 0"))
countMat_d0 <- GetAssayData(cell_subset_d0, slot = "counts")

#get common sig genes between D1vsD0 and D28vsD0
ly6c <- read.csv("/scratch/labs/khatrilab/aleewq/sig_genes_csv/mDc+_day0_vs_28_wilcox_logfc_0.05_cicero_seur_obj_markers.csv")
ly6c_d1 <- read.csv("/scratch/labs/khatrilab/aleewq/sig_genes_csv/mDc+_day0_vs_1_wilcox_logfc_0.05_cicero_seur_obj_markers.csv")

commonGenes <- intersect(ly6c_d1$gene, ly6c$gene)
d1_fc <- ly6c_d1 %>% dplyr::filter(gene %in% commonGenes) %>% dplyr::select(gene,avg_logFC) %>% arrange(desc(avg_logFC))
d28_fc <- ly6c %>% dplyr::filter(gene %in% commonGenes) %>% dplyr::select(gene,avg_logFC)
df <- inner_join(d1_fc,d28_fc, by = "gene")
df <- cbind(df, data.frame("Day0"=rep(0, length(commonGenes))))
colnames(df) <- c("gene", "Day1", "Day28", "Day0")
#plot genes of largest diff between D1vsD0 and D28vsD0
newdf <- df %>% mutate(d1_d28_diff = Day1-Day28) %>% arrange(abs(d1_d28_diff))
newdf_subset <- newdf[1:10,]  #pick top 10 genes (with smallest diff) to plot
df_forPlot <- newdf_subset %>% gather(day, logFC, Day1:Day0)
ggplot(df_forPlot, aes(y=logFC, x=day, group=gene)) + geom_line(aes(color=gene)) + geom_hline(aes(yintercept = 0))

#perform hierarchical clustering of genes based on how 'open' they are on day 28
hclust_matrix <- specific_countMat[commonGenes,]
gene_dist <- dist(hclust_matrix)
gene_hclust <- hclust(gene_dist, method = "complete")
plot(gene_hclust, labels = FALSE)
gene_cluster <- cutree(gene_hclust, k = 4) %>% 
  # turn the named vector into a tibble
  enframe() %>% 
  # rename some of the columns
  rename(gene = name, cluster = value)

##plot mean cicero gene activity trend on D1 and D28 for genes in each cluster
cts_mean <- apply(specific_countMat[commonGenes,], 1, function(x) sum(x)/length(x))  #get mean expression of each gene day28
cts_mean_d1 <- apply(countMat_d1[commonGenes,], 1, function(x) sum(x)/length(x))  #get mean expression of each gene on day1
cts_mean_d0 <- apply(countMat_d0[commonGenes,], 1, function(x) sum(x)/length(x))  #get mean expression of each gene on day1
cts_mean <- as.data.frame(cts_mean)
genes <- rownames(cts_mean)
cts_mean_df <- cts_mean %>% mutate(gene = genes)
cts_mean_d1 <- as.data.frame(cts_mean_d1)
genes2 <- rownames(cts_mean_d1)
cts_mean_df_d1 <- cts_mean_d1 %>% mutate(gene = genes2)
cts_mean_d0 <- as.data.frame(cts_mean_d0)
genes3 <- rownames(cts_mean_d0)
cts_mean_df_d0 <- cts_mean_d0 %>% mutate(gene = genes3)
cts_mean_d1d28d0 <- Reduce(inner_join, list(cts_mean_df, cts_mean_df_d1, cts_mean_df_d0))
cts_cluster <- inner_join(cts_mean_d1d28d0, gene_cluster, by="gene")
cts_cluster_group <- cts_cluster %>% arrange(cluster)
colnames(cts_cluster_group) <- c("day28", "gene", "day1", "day0" , "cluster")
cluster_forPlot <- cts_cluster_group %>% gather(day, mean_cicero_score, c(day28, day1, day0))
ggplot(cluster_forPlot, aes(day, mean_cicero_score)) +
  geom_line(aes(group = gene)) +
  facet_wrap(~cluster)


subset_countMat <- subset(data, idents = "ly6c+, Day 28")
countMat <- GetAssayData(subset_countMat, slot = "counts")
countMat_geneInterest <- countMat[commonGenes,]
geneSeuratObj <- CreateSeuratObject(counts = t(as.matrix(countMat_geneInterest)))
geneSeuratObj <- FindVariableFeatures(object = geneSeuratObj)
geneSeuratObj <- ScaleData(object = geneSeuratObj)
geneSeuratObj <- RunPCA(object = geneSeuratObj)
geneSeuratObj <- FindNeighbors(object = geneSeuratObj)
geneSeuratObj <- FindClusters(object = geneSeuratObj)
geneSeuratObj <- RunTSNE(object = geneSeuratObj, check_duplicates = FALSE)
DimPlot(object = geneSeuratObj, reduction = "tsne", group.by = 'ident')
DoHeatmap(geneSeuratObj, group.by = 'ident', label=TRUE)
geneCluster_seurat <- as.data.frame(geneSeuratObj$seurat_clusters)
fwrite(geneCluster_seurat, file="/home/aleewq/geneClusters.csv", row.names = T)

colnames(geneCluster_seurat) <- "cluster"
genes4 <- rownames(geneCluster_seurat)
geneCluster_seurat <- geneCluster_seurat %>% mutate(gene = genes4)
# cts_cluster <- inner_join(cts_mean_d1d28d0, geneCluster_seurat, by="gene")
# cts_cluster_group <- cts_cluster %>% arrange(cluster)
# colnames(cts_cluster_group) <- c("day28", "gene", "day1", "day0" , "cluster")
# cluster_forPlot <- cts_cluster_group %>% gather(day, mean_cicero_score, c(day28, day1, day0))
# 
# ggplot(cluster_forPlot, aes(day, mean_cicero_score)) +
#   geom_line(aes(group = gene)) +
#   geom_line(stat = "summary", fun.y = "median", colour = "brown", size = 1.5, aes(group = 1)) +
#   facet_wrap(~cluster)

#plot genes of largest diff between D1vsD0 and D28vsD0
geneInterest <- c("Oasl1", "Ptk2b", "Zbp1", "Apobec3", "Irgm1", "Ifit2", "Ly9", "Tap2", "Oas1f", "Nlrc5", 
                  "Ifit1bl1", "Ifit1", "Il6", "Tlr7", "Aif1", "Cxcl10", "Ptafr", "Mx2", "Tnfrsf1a", "Trafd1", "Themis2", 
                  "Oas1c", "Oasl2", "Mx1", "Ifit1bl2", "Ifit3", "Rsad2", "Slfn8")
newdf <- df %>% mutate(d1_d28_diff = Day1-Day28) %>% inner_join(geneCluster_seurat, by="gene") %>% dplyr::filter(gene%in%geneInterest)
df_forPlot <- newdf %>% gather(day, logFC, Day1:Day0)
ggplot(df_forPlot, aes(y=logFC, x=day)) + geom_line(aes(group=gene)) + geom_hline(aes(yintercept = 0)) +
  geom_line(stat = "summary", fun.y = "median", colour = "brown", size = 1.5, aes(group = 1)) +
  facet_wrap(~cluster)


