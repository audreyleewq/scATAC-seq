library(data.table)
library(Matrix)
library(SummarizedExperiment)
library(parallel)
library(Seurat) 
library(dplyr)
library(ggplot2)
library(gridExtra)
library(gmodels)
set.seed(7)

# 3M052 - full mat no subsetting
cicero_seur_obj <- readRDS("/scratch/labs/khatrilab/aleewq/data/fresh_3m052_ciceroobj_logGA_allcells.RDS")
plotUMAP_allcells <- readRDS("/scratch/labs/khatrilab/aleewq/data/plotUMAP_allcells_forplot_withmdcs_andNKs.Rds")

# YF
cicero_seur_obj <- readRDS("/scratch/labs/khatrilab/aleewq/yellow_fever_data/all_yellow_samps_v2_logGA_in_dat_overwrite_cicero_seur_obj.rds")
plotUMAP_allcells <- readRDS("/scratch/labs/khatrilab/aleewq/data/yellowfever_plotUMAP_allcells_forplot_withmdcs_andNKs.Rds")

color_df <- fread("/scratch/labs/khatrilab/aleewq/csvs/cell_color_converter_v3.csv")
manual_scale_colors_clusts <- setNames(color_df$cell_clusts_colors, color_df$cell_clusts)
manual_scale_colors_cell_times <- setNames(color_df$cell_times_colors, color_df$cell_times)

ident_mat <- data.frame(ids = Idents(cicero_seur_obj),
                        barcode = colnames(cicero_seur_obj@assays$RNA@counts))

ident_mat <- ident_mat %>% left_join(plotUMAP_allcells)

cicero_seur_obj@meta.data$cell[1:5]

ident_mat$justcelltype <- sapply(strsplit(ident_mat$cell_2, ", "), "[", 1)

objmat <- cicero_seur_obj@assays$RNA@counts
dim(objmat)
all(colnames(objmat) == ident_mat$barcode)

inds <- which(ident_mat$justcelltype == "dcs")

obj_mat_dcs <- objmat[,inds]
ident_mat_dcs <- ident_mat[inds,]

####

set.seed(7)
obj_mat_dcs_pcs <- fast.prcomp(obj_mat_dcs)
saveRDS(file = "/labs/khatrilab/scottmk/ataq/data/obj_mat_dcs_pcs_fastprcomp_3m052.RDS", obj_mat_dcs_pcs)
# obj_mat_dcs_pcs <- readRDS("/scratch/labs/khatrilab/aleewq/data/obj_mat_dcs_pcs_fastprcomp_3m052.RDS")  #I dont have this 
the_pcs <- obj_mat_dcs_pcs$rotation

srt_dcs <- CreateSeuratObject(obj_mat_dcs) 
srt_dcs <- ScaleData(srt_dcs)
srt_dcs@assays$RNA@scale.data <- as.matrix(srt_dcs@assays$RNA@counts)  #don't actually want data scaled within DCs but need to scale data to get slot
srt_dcs <- FindVariableFeatures(srt_dcs)
srt_dcs <- RunPCA(srt_dcs, features = rownames(obj_mat_dcs))
# overwrite yucky pcs
srt_dcs@reductions$pca@cell.embeddings <- the_pcs[,1:50]
srt_dcs <- FindNeighbors(srt_dcs, dims = 1:10)
srt_dcs <- FindClusters(srt_dcs, resolution = .15)
srt_dcs <- RunTSNE(srt_dcs, dims = 1:10, perplexity = 100)
TSNEPlot(srt_dcs, pt.size=0.5) + labs(color='dc_clust')
saveRDS(file = "/scratch/labs/khatrilab/aleewq/yellow_fever_data/YF_srt_dc_subclusters.RDS", srt_dcs)

# TSNE plot by day
srt_dcs@meta.data$day <- ident_mat_dcs$day
Idents(srt_dcs) <- srt_dcs@meta.data$day
TSNEPlot(srt_dcs, pt.size=0.1) + labs(color='Day')

Idents(srt_dcs) <- srt_dcs@meta.data$RNA_snn_res.0.15 
louvain_ids <- data.frame(dc_subclust = Idents(srt_dcs),
                          barcode = names(Idents(srt_dcs)))

all_dc_marks <- FindAllMarkers(srt_dcs)
fwrite(file = "/scratch/labs/khatrilab/aleewq/yellow_fever_csvs/YF_dc_subclusters_markers.csv", all_dc_marks)

FeaturePlot(srt_dcs, features = "Sirpa") 

#Put labels on Tsne plot
levels(srt_dcs)
new_cluster_ids <- c("CD8- cDC", "migratory DCs", "CD8+/CD103+ cDC", "CD8- cDC", "CD8- cDC", "CD8- cDC")
names(new_cluster_ids) <- levels(srt_dcs)
srt_dcs_newIdents <- RenameIdents(srt_dcs, new_cluster_ids)
TSNEPlot(srt_dcs_newIdents, label=TRUE, pt.size=0.5) 

############## PLOTTING ################
####Plot heatmap of different clusters
library(ComplexHeatmap)

all_dc_marks <- all_dc_marks[-grep("Rik", all_dc_marks$gene),]
all_dc_marks$pct.1_minus_pct.2 <- all_dc_marks$pct.1 - all_dc_marks$pct.2
all_dc_marks <- all_dc_marks[order(all_dc_marks$pct.1_minus_pct.2, decreasing = TRUE),]  #arrange top genes based on pct
all_dc_marks_top <- rbindlist(lapply(split(as.data.table(all_dc_marks), by = "cluster"), function(x) x[1:10,]))
all_dc_marks_top$cluster <- as.character(all_dc_marks_top$cluster)
all_dc_marks_top <- all_dc_marks_top[order(all_dc_marks_top$cluster),]  #order cluster
genes_top <- all_dc_marks_top$gene
genes_ind <- which(rownames(obj_mat_dcs) %in% genes_top)
hmplot <- obj_mat_dcs[genes_ind,]

#set cell barcodes as clusters
plotUMAP_dcs <- as.data.frame(Embeddings(srt_dcs, reduction = "tsne"))
plotUMAP_dcs$barcode <- colnames(srt_dcs@assays$RNA@counts)
plotUMAP_dcs$cluster <- srt_dcs@meta.data$RNA_snn_res.0.15
plotUMAP_dcs_orded <- plotUMAP_dcs[order(plotUMAP_dcs$cluster),]
test <- plotUMAP_allcells %>% dplyr::select(barcode, day, cell_2)
plotUMAP_dcs <- plotUMAP_dcs %>% inner_join(test, by="barcode")
saveRDS(file = "/scratch/labs/khatrilab/aleewq/yellow_fever_data/YF_srt_dcs_plotUMAP.RDS", plotUMAP_dcs)

hmplot <- hmplot[,plotUMAP_dcs_orded$barcode]  #order heatmap matrix according to barcode based on ordered clusters
hmplot <- hmplot[genes_top,]  #order genes based on ordered clusters

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

cols = gg_color_hue(10)

ha = HeatmapAnnotation(clust = plotUMAP_dcs_orded$cluster,
                       col = list(clust = c("0" = cols[1],
                                            "1" = cols[2],
                                            "2" = cols[5],
                                            "3" = cols[6],
                                            "4" = cols[7],
                                            "5" = cols[8])))

pdf("/scratch/labs/khatrilab/aleewq/plots/YF_dcs_subclusters_heatmap.pdf", width = 8, height = 6)
Heatmap(as.matrix(hmplot), show_row_dend = FALSE, show_column_dend = FALSE, show_column_names = FALSE,
        cluster_rows = FALSE, cluster_columns = FALSE,  top_annotation = ha, row_names_gp = gpar(fontsize=7))
dev.off()

######### D28vsD0 differences within cluster #########
identical(rownames(srt_dcs@meta.data), plotUMAP_dcs$barcode) #check samples are in the same order
srt_dcs@meta.data$day <- plotUMAP_dcs$day
saveRDS(file = "/scratch/labs/khatrilab/aleewq/data/3M052_srt_dcs_withclusts.RDS", srt_dcs)

dc_subset = "0"  #dc clust
dc_subset_seurat <- subset(srt_dcs, idents=dc_subset)  #migratory DCs
length(which(dc_subset_seurat$day == "Day 28"))
length(which(dc_subset_seurat$day == "Day 0"))
length(which(dc_subset_seurat$day == "Day 1"))

#Find differences at D28vsD0
Idents(dc_subset_seurat) <- "day"
LOG_FC_THRESH = 0.5
TEST_USE = "wilcox"
day28_day0_markers <- FindMarkers(dc_subset_seurat, ident.1 = "Day 28", ident.2 = "Day 0",
                                    logfc.threshold = LOG_FC_THRESH, test.use = TEST_USE, min.pct = 0, min.diff.pct = 0,
                                    min.cells.feature = 2) # , logfc.threshold = 0.05
day28_day0_markers$pct.1_minus_pct.2 <- day28_day0_markers$pct.1 - day28_day0_markers$pct.2
day28_day0_markers <- day28_day0_markers[order(abs(day28_day0_markers$pct.1_minus_pct.2), decreasing = TRUE),]
day28_day0_markers$gene <- rownames(day28_day0_markers)
day28_day0_markers_sig <- day28_day0_markers[which(day28_day0_markers$p_val_adj < 0.05),]
fwrite(day28_day0_markers_sig, paste0("/home/aleewq/", "dc_clust", dc_subset, "_day28_vs_0_", TEST_USE, "_logfc_", LOG_FC_THRESH, "_srt_dcs_markers.csv"))

#Find differences at D1vsD0
day1_day0_markers <- FindMarkers(dc_subset_seurat, ident.1 = "Day 1", ident.2 = "Day 0",
                                  logfc.threshold = LOG_FC_THRESH, test.use = TEST_USE, min.pct = 0, min.diff.pct = 0,
                                  min.cells.feature = 2) # , logfc.threshold = 0.05
day1_day0_markers$pct.1_minus_pct.2 <- day1_day0_markers$pct.1 - day1_day0_markers$pct.2
day1_day0_markers <- day1_day0_markers[order(abs(day1_day0_markers$pct.1_minus_pct.2), decreasing = TRUE),]
day1_day0_markers$gene <- rownames(day1_day0_markers)
day1_day0_markers_sig <- day1_day0_markers[which(day1_day0_markers$p_val_adj < 0.05),]
fwrite(day1_day0_markers_sig, paste0("/home/aleewq/", "dc_clust", dc_subset, "_day1_vs_0_", TEST_USE, "_logfc_", LOG_FC_THRESH, "_srt_dcs_markers.csv"))


## Look at boxplot for gene
gene_nm <- "Ifit3b"
dim(obj_mat_dcs)

plotUMAP_dcs$geneInterest <- obj_mat_dcs[which(rownames(obj_mat_dcs) == gene_nm),]
plotUMAP_nozero <- plotUMAP_dcs[which(plotUMAP_dcs$geneInterest > 0),]

pzeros <- ggplot(plotUMAP_dcs, aes(day, geneInterest))  + 
  geom_boxplot(color = "red")+ geom_jitter(height = 0, alpha = .5) + 
  ylab(gene_nm) + ggtitle("With Zeros - DCs 3M")
pnozeros <- ggplot(plotUMAP_nozero, aes(day, geneInterest)) + 
  geom_boxplot(color = "red") + geom_jitter(height = 0, alpha = .5) +
  ylab(gene_nm)+ ggtitle("Without Zeros - DCs 3M")

pdf(paste("/home/aleewq/", "dcs_", gene_nm, ".pdf"), width = 8, height = 6)
grid.arrange(pzeros, pnozeros)
dev.off()


####################
#
# find markers on obj
#
####################

#clust1_umap <- uwot::umap(t(as.matrix(clust1)))
set.seed(7)

obj_mat_dcs_pcs <- fast.prcomp(obj_mat_dcs)

saveRDS(file = "/labs/khatrilab/scottmk/ataq/data/obj_mat_dcs_pcs_fastprcomp_3m052.RDS", obj_mat_dcs_pcs)

obj_mat_dcs_pcs <- readRDS("/labs/khatrilab/scottmk/ataq/data/obj_mat_dcs_pcs_fastprcomp_3m052.RDS")
the_pcs <- obj_mat_dcs_pcs$rotation
dim(the_pcs)
the_pcs[1:5,1:5]
obj_mat_dcs[1:5,1:5]

perplex <- as.list(c(5, 15, 30, 50, 100, 300))
sub_tsne_l <- lapply(perplex , function(x){
  res = Rtsne::Rtsne(the_pcs[,1:10], pca = F, perplexity = x)
  return(res)})

#saveRDS(file = "/labs/khatrilab/scottmk/ataq/data/tsne_3m052_dcs_perplex_l.RDS", sub_tsne_l)

sub_tsne_l <- readRDS("/labs/khatrilab/scottmk/ataq/data/tsne_3m052_dcs_perplex_l.RDS")
sub_tsne <- sub_tsne_l[[5]]
#sub_tsne <- Rtsne::Rtsne(t(as.matrix(the_pcs)))
#saveRDS(file = "/labs/khatrilab/scottmk/ataq/data/tsne_3m052_dcs.RDS", sub_tsne)
sub_tsne <- sub_tsne$Y
ident_mat_dcs$tSNE_1_dcs <- sub_tsne[,1]
ident_mat_dcs$tSNE_2_dcs <- sub_tsne[,2]

###################################################
# other clusters tried and failed
###################################################
#### kmeans
# obj_mat_dcs_pcs <- readRDS( "/labs/khatrilab/scottmk/ataq/data/obj_mat_dcs_pcs_fastprcomp_3m052.RDS")
# rot <- obj_mat_dcs_pcs$rotation
# dim(rot)
# rot[1:5,1:5]
# rot_10 <- rot[,1:10]
# set.seed(7)
# irisCluster <- kmeans(rot_10, 6, nstart = 1000)
# ident_mat_dcs$kmeans_clust <- irisCluster$cluster
# 
# ### heirarch
# correlationmatrix = cor(t(ident_mat_dcs[,c("tSNE_1_dcs", "tSNE_2_dcs")]))
# distancematrix <- cor2dist(correlationmatrix)
# clusters <- hclust(as.dist(distancematrix))
# clusterCut <- cutree(clusters, 6)
# ident_mat_dcs$hier_clust <- clusterCut


p0 <- ggplot(ident_mat_dcs[which(ident_mat_dcs$day == "Day 0"),], aes(tSNE_1,tSNE_2,color=cell_2)) + geom_point(size = 1) + 
  theme_bw() + xlab("tSNE1") + ylab("tSNE2")+  guides(colour = guide_legend(override.aes = list(size=3))) +
  theme(text = element_text(size=20)) +  scale_colour_manual(values = manual_scale_colors_cell_times)

p28 <- ggplot(ident_mat_dcs[which(ident_mat_dcs$day == "Day 28"),], aes(tSNE_1,tSNE_2,color=cell_2)) + geom_point(size = 1) + 
  theme_bw() + xlab("tSNE1") + ylab("tSNE2")+  guides(colour = guide_legend(override.aes = list(size=3))) +
  theme(text = element_text(size=20)) +  scale_colour_manual(values = manual_scale_colors_cell_times)

grid.arrange(p0, p28)

p1 <- ggplot(ident_mat_dcs, aes(tSNE_1_dcs,tSNE_2_dcs,color=cell_2)) + geom_point(size = 1) + 
  theme_bw() + xlab("tSNE1") + ylab("tSNE2")+  guides(colour = guide_legend(override.aes = list(size=3))) +
  theme(text = element_text(size=20)) +  scale_colour_manual(values = manual_scale_colors_cell_times)

p2 <- ggplot(ident_mat_dcs, aes(tSNE_1_dcs,tSNE_2_dcs,color=factor(kmeans_clust))) + geom_point(size = 1) + 
  theme_bw() + xlab("tSNE1") + ylab("tSNE2")+  guides(colour = guide_legend(override.aes = list(size=3))) +
  theme(text = element_text(size=20)) 

p2 <- ggplot(ident_mat_dcs, aes(tSNE_1_dcs,tSNE_2_dcs,color=factor(kmeans_clust))) + geom_point(size = 1) + 
  theme_bw() + xlab("tSNE1") + ylab("tSNE2")+  guides(colour = guide_legend(override.aes = list(size=3))) +
  theme(text = element_text(size=20)) 


p2
pdf("/labs/khatrilab/scottmk/ataq/images/subcluster_dcs_v2.pdf")
grid.arrange(p1,p3)
#TSNEPlot(srt_dcs, pt.size = 1)
dev.off()

p2 <- ggplot(ident_mat_dcs, aes(tSNE_1,tSNE_2,color=cell_2)) + geom_point(size = .75) + 
  theme_bw() + xlab("tSNE1") + ylab("tSNE2")+  guides(colour = guide_legend(override.aes = list(size=3))) +
  theme(text = element_text(size=20)) +  scale_colour_manual(values = manual_scale_colors_cell_times)
p2

ident_mat_more <- ident_mat_dcs %>% 
  left_join(louvain_ids)

p3 <- ggplot(ident_mat_more, aes(tSNE_1_dcs,tSNE_2_dcs,color=dc_subclust)) + geom_point(size = 1) + 
  theme_bw() + xlab("tSNE1") + ylab("tSNE2")+  guides(colour = guide_legend(override.aes = list(size=3))) +
  theme(text = element_text(size=20)) 
p3

#################

# raw subcluster
# to subcluster raw atac you have to use seurat lsi or it breaks
obj <- readRDS("/labs/khatrilab/scottmk/ataq/data/ALLSAMPS_scATAC-Summarized-Experiment.rds")
obj_mat <- obj@assays$data@listData$counts
obj_mat_dcs <- obj_mat[,which(colnames(obj_mat) %in% ident_mat_dcs$barcode)]
dim(obj_mat_dcs)
rownames(obj_mat_dcs) <- seq(1:nrow(obj_mat_dcs))
obj_mat_dcs[1:5,1:5]

raw_dc_seurat <- CreateSeuratObject(obj_mat_dcs)
raw_dc_seurat = FindVariableFeatures(raw_dc_seurat, do.plot = F, display.progress = FALSE)
raw_dc_seurat@assays$RNA@scale.data <- raw_dc_seurat@assays$RNA@data

raw_dc_seurat <- RunPCA(object = raw_dc_seurat, verbose = FALSE, )
cicero_seur_obj <- FindNeighbors(object = cicero_seur_obj, dims = 1:7)
cicero_seur_obj <- FindClusters(object = cicero_seur_obj, resolution = 0.2)
cicero_seur_obj <- RunTSNE(object = cicero_seur_obj, dims = 1:7, check_duplicates = FALSE)

