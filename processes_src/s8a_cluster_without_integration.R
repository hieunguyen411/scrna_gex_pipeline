s8a.cluster.wo.integration <- function(s.obj, 
                                      path.to.output, 
                                      save.RDS.s8a,
                                      PROJECT,
                                      num.PCA,
                                      num.PC.used.in.UMAP,
                                      num.PC.used.in.Clustering,
                                      cluster.resolution = 0.5,
                                      my_random_seed = 42,
                                      umap.method = "uwot"){
  chosen.assay <- "RNA"
  DefaultAssay(s.obj) <- chosen.assay
  
  s.obj <- RunPCA(s.obj, npcs = num.PCA, verbose = FALSE, reduction.name=sprintf("%s_PCA", chosen.assay))
  
  s.obj <- RunUMAP(s.obj, reduction = sprintf("%s_PCA", chosen.assay), 
                   dims = 1:num.PC.used.in.UMAP, reduction.name=sprintf("%s_UMAP", chosen.assay),
                   seed.use = my_random_seed, umap.method = "uwot")
  
  # clustering 
  s.obj <- FindNeighbors(s.obj, reduction = sprintf("%s_PCA", chosen.assay), dims = 1:num.PC.used.in.Clustering)
  
  s.obj <- FindClusters(s.obj, 
                        resolution = cluster.resolution, random.seed = 0)
  
  if (save.RDS.s8a == TRUE){
    dir.create(file.path(path.to.output, "s8a_output"), showWarnings = FALSE)
    saveRDS(object = s.obj, 
            file = file.path(path.to.output, "s8a_output", 
                             paste0(PROJECT, ".output.s8a.rds")))
  }
  
  return(s.obj)
}