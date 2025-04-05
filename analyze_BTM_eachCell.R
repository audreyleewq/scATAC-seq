library(Seurat)
library(dplyr)
library(tidyverse)
library(biomaRt)
library(org.Hs.eg.db)
library(corrplot)
library(data.table)

#Collapse to BTM level for each cell
data <- readRDS("/scratch/labs/khatrilab/aleewq/raw_data/samps_10-12/samps_10_12mdcs_beforefindallmarkersbyday_cicero_seur_obj.rds")
plotUMAP_allcells_forplot_withmdcs_andNKs <- readRDS("/scratch/labs/khatrilab/aleewq/data/plotUMAP_allcells_forplot_withmdcs_andNKs.Rds")
countMatrix <- GetAssayData(data, slot = "counts")  #get count matrix
colorLS=colorRampPalette(colors = c("blue", "white", "red"))(n = 100)
cell_interest <- "ly6c+, Day 28"
# "NK, Day 0"            "ly6c+, Day 0"         "NK, Day 1"            "ly6c+, Day 1"         "NK, Day 28"           "ly6c+, Day 28"       
# "cd11b+pdca1+, Day 1"  "cd11b+pdca1+, Day 28" "pdcs, Day 0"          "pdcs, Day 1"          "pdcs, Day 28"         "dcs, Day 0"          
# "dcs, Day 1"           "dcs, Day 28" 

#Convert genes to entrez ID - first convert mouse genes to human genes then convert human genes to entrezID
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")   
mouse <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
genes_human <- getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = rownames(countMatrix), mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
entrez_gene_lookup=select(org.Hs.eg.db, genes_human$HGNC.symbol, c("ENTREZID"), "ALIAS")

#Remove genes with multiple entrez IDs
entrez_gene_lookup=entrez_gene_lookup[!(duplicated(entrez_gene_lookup[,2]) | duplicated(entrez_gene_lookup[,2], fromLast=TRUE)),]
#Remove entrez IDs with multiple genes
entrez_gene_lookup=entrez_gene_lookup[!(duplicated(entrez_gene_lookup[,1]) | duplicated(entrez_gene_lookup[,1], fromLast=TRUE)),]
#generate eset with human genes instead of mouse genes
mouse_human_genes <- genes_human %>% dplyr::filter(HGNC.symbol %in% entrez_gene_lookup[,1])  #remove duplicated human genes
eset <- as.data.frame(countMatrix) %>% dplyr::mutate(MGI.symbol = rownames(countMatrix)) %>% inner_join(mouse_human_genes, by = "MGI.symbol")
rownames(eset) <- eset$HGNC.symbol  
eset <- eset %>% dplyr::select(-c(MGI.symbol, HGNC.symbol))
entrez_list=as.character(entrez_gene_lookup[match(rownames(eset),entrez_gene_lookup[,1]),2])  #get entrezID of genes that are in eset
#Load BTMs
load('/scratch/labs/khatrilab/aleewq/BTM_for_GSEA_20131008_geneID.RData')
#Remove BTMs which have no matching genes in dataset
BTM_list=BTM_list[!sapply(BTM_list, function(x) is_empty(intersect(x,entrez_list)))]
#Collapse - match entrezID in BTM_list to entrezID of genes (arrange according to rownames of eset) - find arithmetic mean for each sample
exp_BTM=do.call(rbind, lapply(BTM_list, function(x) colMeans(eset[na.omit(match(x,entrez_list)),],na.rm=TRUE)))
BTM_interest <- c("M75 - antiviral IFN signature", "M150 - innate antiviral response", "M86.0 - chemokines and inflammatory molecules in myeloid cells", 
                  "M86.1 - proinflammatory dendritic cell, myeloid cell response", "M16 - TLR and inflammatory signaling", "M53 - inflammasome receptors and signaling",
                  "M29 - proinflammatory cytokines and chemokines", "M33 - inflammatory response", "M83 - enriched in naive and memory B cells")
df_sampleID <- plotUMAP_allcells_forplot_withmdcs_andNKs %>% dplyr::select(sample, barcode, cell_2, day)
barcodes_interest <- df_sampleID %>% dplyr::filter(cell_2 == cell_interest) %>% dplyr::select(barcode)  #select cells you're interested to look at
exp_BTM_subset <- t(exp_BTM[BTM_interest, barcodes_interest[,1]])
cor_matrix <- cor(exp_BTM_subset, method = "pearson", use="complete.obs")
corrplot(cor_matrix, method="circle", col=colorLS, tl.cex = 0.6, tl.srt = 45, tl.col='black', win.asp = 1)

#check genes in BTM modules - creates a matrix of number of overlapping genes between the BTMs
BTM_list_subset <- BTM_list[BTM_interest]
BTM_similarity <- df <- data.frame(matrix(ncol = length(BTM_list_subset), nrow = length(BTM_list_subset)))
for(i in 1:length(BTM_interest)) {
  similar <- sapply(1:length(BTM_list_subset), function(x) length(intersect(BTM_list_subset[[x]], BTM_list_subset[[i]])))
  BTM_similarity[i,] <- similar
}
colnames(BTM_similarity) <- BTM_interest
rownames(BTM_similarity) <- BTM_interest
fwrite(BTM_similarity, file = "BTM_interest_similarity.csv", row.names = T)

########

##split BTM modules and do correlation
common_antiviral_genes <- intersect(BTM_list_subset[[1]], BTM_list_subset[[2]])
unique_antiviral_genes_1 <- setdiff(BTM_list_subset[[1]], BTM_list_subset[[2]])
unique_antiviral_genes_2 <- setdiff(BTM_list_subset[[2]], BTM_list_subset[[1]])
BTM_inflammatory <- BTM_interest[-c(1:2, length(BTM_interest))]
inflammatory_BTMs <- BTM_list_subset[BTM_inflammatory]
inflammatory_genes <- c()
for(i in 1:length(inflammatory_BTMs)) {
    genes <- inflammatory_BTMs[[i]]
    inflammatory_genes <- union(inflammatory_genes, genes)
}
B_cell_genes <- BTM_list_subset[[length(BTM_interest)]]
new_BTM_modules <- list("common_antiviral" = common_antiviral_genes, "antiviral_IFN_1" = unique_antiviral_genes_1, 
                        "antiviral_innate_2" = unique_antiviral_genes_2, "inflammatory_response" = inflammatory_genes, "naive_memory_B_cell" = B_cell_genes)
##Collapse genes into newly created BTM modules
exp_BTM=do.call(rbind, lapply(new_BTM_modules, function(x) colMeans(eset[na.omit(match(x,entrez_list)),],na.rm=TRUE)))
df_sampleID <- plotUMAP_allcells_forplot_withmdcs_andNKs %>% dplyr::select(sample, barcode, cell_2, day)
barcodes_interest <- df_sampleID %>% dplyr::filter(cell_2 == cell_interest) %>% dplyr::select(barcode)  #select cells you're interested to look at
exp_BTM_subset <- t(exp_BTM[, barcodes_interest[,1]])
cor_matrix <- cor(exp_BTM_subset, method = "pearson", use="complete.obs")
corrplot(cor_matrix, method="circle", col=colorLS, tl.cex = 0.6, tl.srt = 45, tl.col='black', win.asp = 1)

#check genes in BTM innate antiviral module to see how well correlated they are to each other
test = eset[na.omit(match(BTM_list_subset[[2]], entrez_list)),]
test = t(test)
cor_matrix <- cor(test, method = "pearson", use="complete.obs")
corrplot(cor_matrix, method="circle", col=colorLS, tl.cex = 0.6, tl.srt = 45, tl.col='black', win.asp = 1)

#find BTM genes that are present in eset 
genes_in_eset <- sapply(new_BTM_modules, function(x) rownames(eset)[na.omit(match(x, entrez_list))])

#######

#perform clustering on genes (sig genes)
##perform on ly6c+ cells at D28
library(stats)
library(tibble)
specific_cell_subset <- subset(data, idents = c("ly6c+, Day 28"))
specific_countMat <- GetAssayData(specific_cell_subset, slot = "counts")
sig_genes <- read.csv("/home/aleewq/res_outputsamps_10_12mdcs_ly6c+, Day 28_test_t_ONLY_0_28_logfc_thresh_0.05.csv", header = T)
sig_genes <- sig_genes %>% select(gene)
sig_genes <- as.character(sig_genes[,1])
hclust_matrix <- specific_countMat[sig_genes,]
gene_dist <- dist(hclust_matrix)
gene_hclust <- hclust(gene_dist, method = "complete")
gene_cluster <- cutree(gene_hclust, k = 5) %>% 
  # turn the named vector into a tibble
  enframe() %>% 
  # rename some of the columns
  rename(gene = name, cluster = value)

cts_mean <- apply(specific_countMat[sig_genes,], 1, function(x) sum(x)/length(x))  #get mean expression of each gene
cts_mean <- as.data.frame(cts_mean)
cts_cluster <- cbind(cts_mean, gene_cluster)
cts_cluster_group <- cts_cluster %>% arrange(cluster)
fwrite(cts_cluster_group, file = "samps10_12_cicero_ly6c_D28_hclust_gene_clusters.csv")


