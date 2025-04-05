#saveRDS(file = "/scratch/labs/khatrilab/aleewq/ArchR_dat/Archr_genemat_magicimputed.rds", magicmat)
#saveRDS(file = "/scratch/labs/khatrilab/aleewq/ArchR_dat/cell_info_v2.RDS", df)
library(ggrastr)
library(dplyr)
library(ggplot2)
library(data.table)
library(ArchR)
library(viridis)
library(purrr)
library(cowplot)
df <- readRDS("/scratch/labs/khatrilab/aleewq/ArchR_dat/cell_info_v2.RDS")
magicmat <- readRDS("/scratch/labs/khatrilab/aleewq/ArchR_dat/Archr_genemat_magicimputed.rds")
magicmat_df <- as.data.frame(t(magicmat))
magicmat_df$UMAP1 <- df$UMAP1
magicmat_df$UMAP2 <- df$UMAP2
magicmat_df$Clusters <- df$Clusters
magicmat_df$sampe_bc <- df$sampe_bc
magicmat_df_dcs <- magicmat_df[which(magicmat_df$Clusters %in% c("C12", "C13", "C14", "C15", "C18", "C19", "C1")),]
df_sub_sub <- subset(magicmat_df_dcs, select = c("sampe_bc", "Clusters",  "UMAP1", "UMAP2",
                                                 "Tgm2", "Kremen1", "Ankrd22", "Klf4"
))

df_sub_melt <- melt(df_sub_sub, id.vars = c("sampe_bc", "Clusters", "UMAP1", "UMAP2"))
# this quantile cuts to limit the effects of outliers on the color scheme
df_sub_melt_l <- split(as.data.table(df_sub_melt), by = "variable")
df_sub_melt_l <- lapply(df_sub_melt_l, function(x){
  varx <- x$value
  newx <- ArchR:::.quantileCut(varx)
  x$value <- newx
  return(x)
})
df_sub_melt <- rbindlist(df_sub_melt_l)
df_sub_melt %>% 
  group_split(variable) %>% 
  map(
    ~ggplot(., aes(UMAP1, UMAP2, color = value)) + 
      geom_point_rast(size = .5) +
      scale_color_viridis(option = "plasma" #, 
                          # limits = c(quantile(value, probs = .025),
                          #        quantile(value, probs = .975)), 
                          # oob = scales::squish
      ) + theme_cowplot() +
      theme(legend.position="right") + xlim(-2,13) + ylim(-8,3) +
      theme(legend.title = element_blank()) +
      facet_grid(~ variable, labeller = function(x) label_value(x, multi_line = FALSE))
  ) %>% 
  plot_grid(plotlist = ., align = 'hv', ncol = 2)


## Look at cell info
df_subset <- df[which(df$Clusters %in% c("C12", "C13", "C14", "C15", "C18", "C19", "C1")),]

ggplot(df_subset, aes(x=UMAP1, y=UMAP2, color=source)) + 
  geom_point_rast(size = .1) + theme_cowplot() + 
  theme(legend.position="right") + xlim(-2,13) + ylim(-8,3) +
  theme(legend.title = element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=2))) + scale_color_manual(values = c("3M"="red", "YF"="gold"))

ggplot(df_subset, aes(x=UMAP1, y=UMAP2, color=day)) + 
  geom_point_rast(size = .1) + theme_cowplot() + 
  theme(legend.position="right") + xlim(-2,13) + ylim(-8,3) +
  theme(legend.title = element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=2))) + 
  scale_color_manual(values = c("Day 0"="mediumorchid4", "Day 1"="chocolate1", "Day 28"="lightskyblue")) 
