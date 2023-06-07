s8a.cluster.wo.integration <- function(s.obj, 
                                      path.to.output, 
                                      save.RDS.s8a,
                                      PROJECT,
                                      num.PCA,
                                      num.PC.used.in.UMAP,
                                      num.PC.used.in.Clustering,
                                      cluster.resolution = 0.5,
                                      my_random_seed = 42,
                                      umap.method = "uwot",
                                      genes.to.not.run.PCA = NULL,
                                      pca_reduction_name = NULL,
                                      umap_reduction_name = NULL){
  chosen.assay <- "RNA"
  DefaultAssay(s.obj) <- chosen.assay
  
  if (is.null(genes.to.not.run.PCA) == TRUE){
    s.obj <- RunPCA(s.obj, npcs = num.PCA, verbose = FALSE, reduction.name=sprintf("%s_PCA", chosen.assay))
    s.obj <- RunUMAP(s.obj, reduction = sprintf("%s_PCA", chosen.assay), 
                     dims = 1:num.PC.used.in.UMAP, reduction.name=sprintf("%s_UMAP", chosen.assay),
                     seed.use = my_random_seed, umap.method = "uwot")
    # clustering 
    s.obj <- FindNeighbors(s.obj, reduction = sprintf("%s_PCA", chosen.assay), dims = 1:num.PC.used.in.Clustering)
    s.obj <- FindClusters(s.obj, resolution = cluster.resolution, random.seed = 0)
    
  } else {
    if(is.null(pca_reduction_name) == TRUE | is.null(umap_reduction_name) == TRUE){
      stop("When genes.to.not.run.PCA is NULL, pca_reduction_name and umap_reduction_name must be used!")
    } else {
      pca.genes <- setdiff(VariableFeatures(s.obj), genes.to.not.run.PCA)
      print("####################################################################")
      print(sprintf("Running PCA with %s genes, after removing %s genes", length(pca.genes), length(genes.to.not.run.PCA)))
      print("####################################################################")
      s.obj <- RunPCA(s.obj, npcs = num.PCA, verbose = FALSE, reduction.name = pca_reduction_name, 
                      features = pca.genes)
      s.obj <- RunUMAP(s.obj, reduction = pca_reduction_name, 
                       dims = 1:num.PC.used.in.UMAP, reduction.name = umap_reduction_name,
                       seed.use = my_random_seed, umap.method = "uwot")
      # clustering 
      s.obj <- FindNeighbors(s.obj, reduction = pca_reduction_name, dims = 1:num.PC.used.in.Clustering)
      s.obj <- FindClusters(s.obj, resolution = cluster.resolution, random.seed = 0)
    }
  }
  
  if (save.RDS.s8a == TRUE){
    dir.create(file.path(path.to.output, "s8a_output"), showWarnings = FALSE)
    saveRDS(object = s.obj, 
            file = file.path(path.to.output, "s8a_output", 
                             paste0(PROJECT, ".output.s8a.rds")))
  }
  
  return(s.obj)
}