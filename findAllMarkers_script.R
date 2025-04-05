library(Seurat)
library(dplyr)
library(fgsea)

# Differential expression
master_file <- "samps_10_12mdcs"
data <- readRDS("/scratch/labs/khatrilab/aleewq/raw_data/samps_10-12/samps_10_12mdcs_beforefindallmarkersbyday_cicero_seur_obj.rds")
plotUMAP_allcells_forplot_withmdcs_andNKs <- readRDS("/scratch/labs/khatrilab/aleewq/data/plotUMAP_allcells_forplot_withmdcs_andNKs.Rds")

#Filter cells to include only ly6c+ cells (excluding NK cells) in this case
ly6c_all <- plotUMAP_allcells_forplot_withmdcs_andNKs %>% dplyr::filter(mDc_or_nk == "mDc") %>% dplyr::select(cell_2, big_clust_ids)

#Set Idents to the group you want 
data@meta.data$cell <- ly6c_all$cell_2
data@meta.data$big_clust_ids <- ly6c_all$big_clust_ids
#find differential markers 
LOG_FC_THRESH = 0.05
TEST_USE = "t"
min.pct = -Inf
Idents(data) <- "cell"  #set active ident 
ident_1 = "ly6c+, Day 28"
ident_2 = "ly6c+, Day 0"
day0_day28_markers <- FindMarkers(data, ident.1 = ident_1, ident.2 = ident_2,
                                  logfc.threshold = LOG_FC_THRESH, test.use = TEST_USE,
                                  min.pct = min.pct, return.thresh = 0,
                                  min.cells.feature = 1) 
day0_day28_markers$pct.1_minus_pct.2 <- day0_day28_markers$pct.1 - day0_day28_markers$pct.2
day0_day28_markers <- day0_day28_markers[order(abs(day0_day28_markers$pct.1_minus_pct.2), decreasing = TRUE),]
day0_day28_markers$gene <- rownames(day0_day28_markers)
day0_day28_markers_sig <- day0_day28_markers[which(day0_day28_markers$p_val_adj < 0.05),]
colnames(day0_day28_markers_sig)[3:4] <- c("pct.1 (Day 28)", "pct.2 (Day 0)")
fwrite(day0_day28_markers_sig, paste0("res_output", master_file, "_", ident_1, "_test_", TEST_USE,  "_ONLY_0_28_logfc_thresh_",LOG_FC_THRESH, ".csv"))

#Do GSEA on sig. genes - BTM
d28vsd0_for_gsea <- day0_day28_markers_sig %>% mutate(sign_pval = -log10(p_val_adj)*sign(avg_logFC)) %>% dplyr::select(gene, sign_pval)
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
sig_genes_human <- getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = d28vsd0_for_gsea$gene, mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
genes_not_found <- setdiff(d28vsd0_for_gsea$gene,sig_genes_human$MGI.symbol)  #mouse genes that are excluded during the mouse-human conversion

#join pval to human genes and order by descreasing pval
d28vsd0_human_gsea <- left_join(sig_genes_human, d28vsd0_for_gsea, by= c("MGI.symbol" = "gene"))  #get pval and converted mouse to human geens
d28vsd0_human_gsea <- d28vsd0_human_gsea %>% arrange(desc(sign_pval))
final_d28vsd0_human <- d28vsd0_human_gsea[!duplicated(d28vsd0_human_gsea$HGNC.symbol),]  ##remove duplicated genes

#get pvals and attribute gene names to pvals
ranks = final_d28vsd0_human$sign_pval
attr(ranks,"names")=as.character(final_d28vsd0_human$HGNC.symbol)

#load BTM
load("/home/aleewq/bloodtranscriptionalmodules.RData")
#Perform GSEA on genes ranked by pvals*sign(logFC)
fgsea_result = fgsea(modGenes, ranks, nperm=1000, maxSize=500)  
test <- fgsea_result %>% filter(padj < 0.05) %>% dplyr::select(pathway, padj, NES, leadingEdge)
fwrite(file = "/home/aleewq/gsea_result", test)

GSEA_NES=as.numeric(fgsea_result$NES)
GSEA_p=fgsea_result$padj
p_cutoff=0.05
ind=which(GSEA_p<p_cutoff)
ind_sort=order(GSEA_NES[ind])
GSEA_sig = data.frame(pathway=fgsea_result$pathway[ind[ind_sort]], 
                      NES=GSEA_NES[ind[ind_sort]], 
                      padj = GSEA_p[ind[ind_sort]])
head(GSEA_sig)