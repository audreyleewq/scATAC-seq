#Chromvar
library(chromVAR)
library(dplyr)
library(SummarizedExperiment)
library(motifmatchr)
library(BSgenome.Mmusculus.UCSC.mm10)
library(Seurat)
library(data.table)

#3M052
data <- readRDS("/scratch/labs/khatrilab/aleewq/data/all_samps_chromVAR-Summarized-Experiment_v2.rds")
plotUMAP <- readRDS("/scratch/labs/khatrilab/aleewq/data/plotUMAP_allcells_forplot_withmdcs_andNKs.Rds")
data@colData$cell_2 <- plotUMAP$cell_2  #colData of summarizedExperiment is like pData of eset (saves the phenoData)
chromvar_seurat_obj <- CreateSeuratObject(data@assays$data$deviations)
Idents(chromvar_seurat_obj) <- plotUMAP$cell_2
saveRDS(chromvar_seurat_obj, file="/scratch/labs/khatrilab/aleewq/data/3M052_chromvar_seur_obj.rds")

#Try plotting Tsne calculated from chromvar deviations
# tsne_results <- deviationsTsne(data, threshold = 1.5, perplexity = 10)
tsne_results <- as.data.frame(colData(data)) %>% dplyr::select(UMAP1, UMAP2)  #use the umap values already calculated from genomic peaks before calling chromvar
tsne_plots <- plotDeviationsTsne(data, tsne_results, 
                                 annotation_name = "Hlf",
                                 sample_column = "cell_2", 
                                 shiny = FALSE)
tsne_plots[[1]]  #gives the whole tsne plot
tsne_plots[[2]]  #gives the featureplot



#For reference - this is how you get TF motifs from chromvar and the deviation for each motif (something like zscore-
#how accessible the peaks is relative to the expectation based on equal chromatin accessbility across the cells. normalized by background peak sets)
#After MACS2 peak calling (output bed files) and combining the bed files across all samples
#you get a count matrix with genomic locations - you save it as a summarized experiment (which is the ALLSAMPS-summarized-experiment)
#GRanges is a vector that stores the genomic location and score - essentially like a way to visualize your bed files
#SummarizedExperiment is something like eset - allows you to store the count matrix and visualize it easily
# se <- readRDS("/scratch/labs/khatrilab/aleewq/data/ALLSAMPS_scATAC-Summarized-Experiment.rds")
# head(assay(se, "counts"))
# genome <- BSgenome.Mmusculus.UCSC.mm10
# se <- addGCBias(se, genome = genome)
# matches <- matchMotifs(mouse_pwms_v1, rowRanges(se), genome = "BSgenome.Mmusculus.UCSC.mm10")
# dev <- computeDeviations(object = se, annotations = matches)

#YF
# data <- readRDS("/scratch/labs/khatrilab/aleewq/yellow_fever_data/all_yellow_samps_v1_chromVAR-Summarized-Experiment.rds")
# yellowfever_plotUMAP_allcells_forplot_withmdcs_andNKs <- readRDS("/scratch/labs/khatrilab/aleewq/data/yellowfever_plotUMAP_allcells_forplot_withmdcs_andNKs.Rds")
# data@colData$cell_2 <- yellowfever_plotUMAP_allcells_forplot_withmdcs_andNKs$cell_2
# chromvar_seurat_obj <- CreateSeuratObject(data@assays$data$deviations)
# Idents(chromvar_seurat_obj) <- yellowfever_plotUMAP_allcells_forplot_withmdcs_andNKs$cell_2

#Find differential deviation between groups
LOG_FC_THRESH = 0.05
TEST_USE = "wilcox" #wilcox
res <- FindMarkers(chromvar_seurat_obj, ident.1 = "NK, Day 1", ident.2 = "NK, Day 0",
            logfc.threshold = LOG_FC_THRESH, test.use = TEST_USE, min.pct = 0, min.diff.pct = 0,
            min.cells.feature = 2) # , logfc.threshold = 0.05

resfilter <- res %>% mutate(motif = rownames(res)) %>% filter(p_val_adj < 0.05)
fwrite(resfilter, file = "/home/aleewq/3M052_chromvar_NK_D1vsD0_wilcox_p05_fc05.csv")


#####using chromvar differential function - basically it is also a wrapper function for t-test/wilcox

test <- differentialDeviationsModified(data, groups="cell_2", alternative = "two.sided", cells = c("ly6c+, Day 0", "ly6c+, Day 28"))
sigMotifs <- test %>% mutate(motif = rownames(test)) %>% dplyr::filter(p_value_adjusted < 0.05)
sigMotifs <- sigMotifs[order(sigMotifs$p_value_adjusted),]
head(sigMotifs)

cells = c("ly6c+, Day 0", "ly6c+, Day 28")
ind <- which(groups %in% cells) #added
groups <- groups[ind]  #added
inputs <- data[,ind]

tsne_results <- deviationsTsne(inputs, threshold = 1.5, perplexity = 10)
tsne_plots <- plotDeviationsTsne(inputs, tsne_results, 
                                 annotation_name = "ENSMUSG00000025498_LINE1499_Irf7_D", 
                                 sample_column = "cell_2", 
                                 shiny = FALSE)
tsne_plots[[1]]
tsne_plots[[2]]


differentialDeviationsModified <- function(object, 
                                   groups,
                                   alternative = c("two.sided", "less",
                                                   "greater"), 
                                   parametric = TRUE,
                                   cells) {
  stopifnot(is(object,"chromVARDeviations"))
  if (length(groups) == 1 && groups %in% colnames(object@colData)) {
    groups <- object@colData[[groups]]
  } else if (length(groups) != ncol(object)) {
    stop("invalid groups input, must be vector of lench ncol(object) or column",
         " name from colData(object)")
  }
  ind <- which(groups %in% cells) #added
  groups <- groups[ind]  #added
  groups <- as.factor(groups)
  
  print(groups)
  
  alternative <- match.arg(alternative)
  inputs <- deviations(object)
  inputs <- inputs[,ind]  ##added
  
  print(dim(inputs))
  
  if (parametric) {
    if (nlevels(groups) == 2) {
      # t-test
      print("doing t-test")
      p_val <- apply(inputs, 1, t_helper, groups, alternative)
    } else {
      # anova
      p_val <- apply(inputs, 1, anova_helper, groups)
    }
  } else {
    if (nlevels(groups) == 2) {
      # wilcoxon
      p_val <- apply(inputs, 1, wilcoxon_helper, groups, alternative)
    } else {
      # kruskal-wallis
      p_val <- apply(inputs, 1, kw_helper, groups)
    }
  }
  
  p_adj <- p.adjust(p_val, method = "BH")
  return(data.frame(p_value = p_val, p_value_adjusted = p_adj))
}


t_helper <- function(x, groups, alternative) {
  splitx <- split(x, groups)
  return(t.test(splitx[[1]],splitx[[2]],
                alternative = alternative, 
                paired = FALSE,
                var.equal = FALSE)$p.value)
}

anova_helper <- function(x, groups) {
  tmpdf <- data.frame(groups = groups, devs = x)
  res <- oneway.test(devs ~ groups, tmpdf, var.equal = FALSE)
  return(res$p.value)
}

kw_helper <- function(x, groups) {
  tmpdf <- data.frame(groups = groups, devs = x)
  res <- kruskal.test(devs ~ groups, tmpdf)
  return(res$p.value)
}

wilcoxon_helper <- function(x, groups, alternative) {
  splitx <- split(x, groups)
  return(wilcox.test(splitx[[1]], splitx[[2]],
                     alternative = alternative, 
                     paired = FALSE)$p.value)
}




