s8.integration.and.clustering <- function(s.obj, 
                           path.to.output, 
                           save.RDS.s8,
                           PROJECT, 
                           num.dim.integration,
                           num.PCA,
                           num.dim.cluster,
                           cluster.resolution = 0.5
                           ){
  data.list <- SplitObject(s.obj, split.by = "name")
  
  data.list <- lapply(X = data.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)})
  
  k.filter <- 200
  
  anchors <- FindIntegrationAnchors(object.list = data.list, dims = 1:num.dim.integration, scale=F,
                                    k.filter = k.filter)## THIS IS CCA DIMENSIONS
  
  s.obj_inte <- IntegrateData(anchorset = anchors, dims = 1:num.dim.integration, k.weight = k.filter) ## THIS IS PCA DIMENSION
  
  ## keep the order of integration obj
  s.obj_inte <- s.obj_inte[, colnames(s.obj)]
  
  s.obj[['integrated']] <- s.obj_inte[['integrated']]
  
  s.obj@commands <- c(s.obj@commands, s.obj_inte@commands)
  
  s.obj@tools <- c(s.obj@tools, s.obj_inte@tools)
  
  DefaultAssay(s.obj) <- "integrated"
  
  s.obj <- ScaleData(s.obj, verbose = FALSE, features = row.names(s.obj))
  
  s.obj <- RunPCA(s.obj, npcs = num.PCA, verbose = FALSE, reduction.name="INTE_PCA")
  
  s.obj <- RunUMAP(s.obj, reduction = "INTE_PCA", dims = 1:num.PCA, reduction.name="INTE_UMAP")
  
  # clustering 
  s.obj <- FindNeighbors(s.obj, reduction = "INTE_PCA", dims = 1:num.dim.cluster)
  
  s.obj <- FindClusters(s.obj, resolution = cluster.resolution)
  
  if (save.RDS.s8 == TRUE){
    dir.create(file.path(path.to.output, "s8_output"), showWarnings = FALSE)
    saveRDS(object = s.obj, 
            file = file.path(path.to.output, "s8_output", 
                             paste0(PROJECT, ".output.s8.rds")))
  }
  
  return(s.obj)
}