## ArchR analysis - re-clustering clusters after Harmony batch correction

subproj_sub_lymph <- proj_peakdat[which(proj_peakdat$Clusters %in% c("C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9")),]
subproj_sub_lymph <- addImputeWeights(subproj_sub_lymph)
subproj_sub_lymph <- addIterativeLSI(
  ArchRProj = subproj_sub_lymph,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI_subpeak_lymph", 
  iterations = 2, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.2), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:30,
  force = T
)
subproj_sub_lymph <- addHarmony(
  ArchRProj = subproj_sub_lymph,
  reducedDims = "IterativeLSI_subpeak_lymph",
  name = "Harmony",
  groupBy = "source",
  force = T
)
subproj_sub_lymph <- addUMAP(
  ArchRProj = subproj_sub_lymph, 
  reducedDims = "Harmony", 
  name = "UMAP_subpeak_Harmony", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine",
  force = T
)
subproj_sub_lymph <- addClusters(
  input = subproj_sub_lymph,
  reducedDims = "Harmony",
  method = "Seurat",
  name = "Clusters_lymph_harmony",
  resolution = 0.4,
  force = T
)