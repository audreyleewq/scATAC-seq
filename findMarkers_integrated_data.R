library(cicero)
library(data.table)
library(Matrix)
library(GenomicRanges)
library(magrittr)
library(SummarizedExperiment)
library(optparse)
library(yaml)
library(Rcpp)
library(parallel)
library(Seurat)
library(dplyr)
set.seed(1)

#"/scratch/labs/khatrilab/aleewq/data/YF_3M_dat/cicero_log_geneactivity_ArchR_overlap.RDS" - seurat object of cicero data
#"/scratch/labs/khatrilab/aleewq/data/YF_3M_dat/cell_info_cicero_ArchR_overlap.RDS" - meta data about the cells in the seurat object
#columns of the gene score matrices from ArchR is not the same as that from the Seurat object

master_file <- "/scratch/labs/khatrilab/aleewq/data/YF_3M_dat/"
setwd(master_file)
# read in cicero data and make seurat object
cicero_logGA <- readRDS("/scratch/labs/khatrilab/aleewq/data/YF_3M_dat/cicero_log_geneactivity_ArchR_overlap.RDS")
data <- cicero_logGA@assays$data@listData$logGA

#read in phenotype info
plotUMAP_allcells <- readRDS("/scratch/labs/khatrilab/aleewq/data/YF_3M_dat/cell_info_cicero_ArchR_overlap.RDS")
#get 3m052 data from integrated data
source_3m <- plotUMAP_allcells %>% filter(source == "3M") 
barcodes_3m <- source_3m$sample_bc
data_3m <- data[,barcodes_3m]
#get YF data from integrated data
source_yf <- plotUMAP_allcells %>% filter(source == "YF") 
barcodes_yf <- source_yf$sample_bc
data_yf <- data[,barcodes_yf]

dataList <- list(data_3m, data_yf) %>% set_names(c("3M", "YF"))
summary(dataList)

#Create separate seurat obj for 3M and YF data
seuratList <- list()
for(i in 1:length(dataList)) {
  data <- dataList[[i]]
  experiment_name <- names(dataList[i])
  cicero_seur_obj <- CreateSeuratObject(data)
  # cicero_seur_obj <- NormalizeData(cicero_seur_obj) # data already logged
  cicero_seur_obj = FindVariableFeatures(cicero_seur_obj, do.plot = F, display.progress = FALSE)
  cicero_seur_obj = ScaleData(cicero_seur_obj, verbose = FALSE)
  
  # this is needed here with the mdcs to get a reasonable tsne
  # cicero_seur_obj@assays$RNA@data <- Matrix(as.matrix(data), sparse = T)
  # cicero_seur_obj@assays$RNA@scale.data <- as.matrix(data)
  
  cicero_seur_obj <- RunPCA(object = cicero_seur_obj, verbose = FALSE, features = rownames(cicero_seur_obj@assays$RNA@data))
  cicero_seur_obj <- FindNeighbors(object = cicero_seur_obj, dims = 1:7)
  cicero_seur_obj <- FindClusters(object = cicero_seur_obj, resolution = 0.2)
  cicero_seur_obj <- RunTSNE(object = cicero_seur_obj, dims = 1:7, check_duplicates = FALSE)
  
  saveRDS(cicero_seur_obj, paste0("/scratch/labs/khatrilab/aleewq/3M_YF_seur_obj_from_intData/", experiment_name, "_cicero_seur_obj.RDS"))
}

######### findMarker ANALYSIS ############
cicero_seur_obj <- readRDS("/scratch/labs/khatrilab/aleewq/3M_YF_seur_obj_from_intData/YF_cicero_seur_obj.RDS")
dim(GetAssayData(cicero_seur_obj))
TSNEPlot(cicero_seur_obj,label = TRUE,label.pt.size = 0.5)
unique(Idents(cicero_seur_obj))
#check barcodes are identical, set idents to be cell types
# identical(source_3m$sample_bc, colnames(GetAssayData(cicero_seur_obj)))
# Idents(cicero_seur_obj) <- source_3m$cell_2

#Get cell type ident into seurat obj
ident_mat <- data.frame(ids = Idents(cicero_seur_obj),
                        sample_bc = colnames(cicero_seur_obj@assays$RNA@counts))
ident_mat <- ident_mat %>% left_join(source_yf, by="sample_bc")
ident_mat$justcelltype <- sapply(strsplit(ident_mat$cell_2, ", "), "[", 1)
cicero_seur_obj@meta.data$cell_2 <- ident_mat$cell_2
cicero_seur_obj@meta.data$justcelltype <- ident_mat$justcelltype

#change ly6c+ to mDc 
cicero_seur_obj@meta.data$justcelltype <-gsub("ly6c+", "mDc", cicero_seur_obj@meta.data$justcelltype)
cicero_seur_obj@meta.data$cell_2 <-gsub("ly6c+", "mDc", cicero_seur_obj@meta.data$cell_2)

#look at overall numbers of each cell type
table(cicero_seur_obj@meta.data$cell_2)
table(cicero_seur_obj@meta.data$justcelltype)

setwd("/scratch/labs/khatrilab/aleewq/3M_YF_seur_obj_from_intData/")
master_file = "/scratch/labs/khatrilab/aleewq/3M_YF_seur_obj_from_intData/YF_results/"

Idents(cicero_seur_obj) <- "cell_2"

## should add in:
LOG_FC_THRESH = 0
TEST_USE = "wilcox" #t

unique(cicero_seur_obj@meta.data$justcelltype)
# cell_type <- "cd11b+pdca1+"


##DAY 28vs0
for(cell_type in unique(cicero_seur_obj@meta.data$justcelltype)){
  print(cell_type)
  cicero_seur_obj_sub <- subset(cicero_seur_obj, subset = justcelltype == cell_type)
  
  
  ident_day0 <- paste0(cell_type, ", Day 0")
  ident_day1 <- paste0(cell_type, ", Day 1")
  ident_day28 <- paste0(cell_type, ", Day 28")
  
  num_day0 <- length(which(Idents(cicero_seur_obj_sub)== ident_day0))
  num_day28 <- length(which(Idents(cicero_seur_obj_sub)== ident_day28))
  if(num_day0 > 1 & num_day28 > 1){
    day0_day28_markers_roc <- FindMarkers(cicero_seur_obj_sub, ident.1 = ident_day28, ident.2 = ident_day0,
                                          logfc.threshold = LOG_FC_THRESH, test.use = TEST_USE, min.pct = 0, min.diff.pct = 0,
                                          min.cells.feature = 2) # , logfc.threshold = 0.05
    
    day0_day28_markers_roc$pct.1_minus_pct.2 <- day0_day28_markers_roc$pct.1 - day0_day28_markers_roc$pct.2
    day0_day28_markers_roc <- day0_day28_markers_roc[order(abs(day0_day28_markers_roc$pct.1_minus_pct.2), decreasing = TRUE),]
    day0_day28_markers_roc$gene <- rownames(day0_day28_markers_roc)
    
    fwrite(day0_day28_markers_roc, paste0(master_file, cell_type, "_day28_vs_0_", TEST_USE, "_logfc_", LOG_FC_THRESH, "_ALL_cicero_seur_obj_markers.csv"))
    
    day0_day28_markers_sig <- day0_day28_markers_roc[which(day0_day28_markers_roc$p_val_adj < 0.05),]
    colnames(day0_day28_markers_sig)[3:4] <- c(paste0("pct.1 (", ident_day28, ")"),   paste0("pct.2 (", ident_day0, ")"))
    fwrite(day0_day28_markers_sig, paste0(master_file, cell_type, "_day28_vs_0_", TEST_USE, "_logfc_", LOG_FC_THRESH, "_pvaladj_0.05_cicero_seur_obj_markers.csv"))
  }
}

### DAY 1
for(cell_type in unique(cicero_seur_obj@meta.data$justcelltype)){
  print(cell_type)
  
  cicero_seur_obj_sub <- subset(cicero_seur_obj, subset = justcelltype == cell_type)
  length(Idents(cicero_seur_obj_sub))
  unique(Idents(cicero_seur_obj_sub))
  
  ident_day0 <- paste0(cell_type, ", Day 0")
  ident_day1 <- paste0(cell_type, ", Day 1")
  
  num_day0 <- length(which(Idents(cicero_seur_obj_sub)== ident_day0))
  num_day1 <- length(which(Idents(cicero_seur_obj_sub)== ident_day1))
  if(num_day0 > 1 & num_day1 > 1){
    day0_day1_markers_roc <- FindMarkers(cicero_seur_obj_sub, ident.1 = ident_day1, ident.2 = ident_day0,
                                         logfc.threshold = LOG_FC_THRESH, test.use = TEST_USE, min.pct = 0, min.diff.pct = 0,
                                         min.cells.feature = 2) # , logfc.threshold = 0.05
    
    day0_day1_markers_roc$pct.1_minus_pct.2 <- day0_day1_markers_roc$pct.1 - day0_day1_markers_roc$pct.2
    day0_day1_markers_roc <- day0_day1_markers_roc[order(abs(day0_day1_markers_roc$pct.1_minus_pct.2), decreasing = TRUE),]
    day0_day1_markers_roc$gene <- rownames(day0_day1_markers_roc)
    
    fwrite(day0_day1_markers_roc, paste0(master_file, cell_type, "_day1_vs_0_", TEST_USE, "_logfc_", LOG_FC_THRESH, "_ALL_cicero_seur_obj_markers.csv"))
    
    day0_day1_markers_sig <- day0_day1_markers_roc[which(day0_day1_markers_roc$p_val_adj < 0.05),]
    colnames(day0_day1_markers_sig)[3:4] <- c(paste0("pct.1 (", ident_day1, ")"),   paste0("pct.2 (", ident_day0, ")"))
    fwrite(day0_day1_markers_sig, paste0(master_file, cell_type, "_day1_vs_0_", TEST_USE, "_logfc_", LOG_FC_THRESH, "_pvaladj_0.05_cicero_seur_obj_markers.csv"))
  }
}