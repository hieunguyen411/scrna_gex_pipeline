s10.PROGENy_pathway_analysis <- function( s.obj, 
                                          path.to.output, 
                                          progeny.params,
                                          cluster.name, 
                                          save.RDS.s10,
                                          PROJECT){
  Idents(s.obj) <- cluster.name
  
  s.obj <- progeny::progeny(s.obj, scale = progeny.params$scale, 
                                    organism = progeny.params$species, 
                                    top = progeny.params$top, 
                                    perm = progeny.params$perm, 
                                    return_assay = progeny.params$return_assay)
  s.obj <- ScaleData(s.obj, features = row.names(s.obj))
  DefaultAssay(s.obj) <-'progeny'
  all.pathways <- rownames(s.obj@assays$progeny)
  
  Idents(s.obj) <- cluster.name
  
  if (is.factor(s.obj@meta.data[, cluster.name])){
    s.obj@meta.data[, cluster.name] <- droplevels(s.obj@meta.data[, cluster.name])
  }else{
    c_names <-unique(s.obj@meta.data[, cluster.name])
    help_sort_func <- ifelse(all.is.numeric(c_names), as.numeric, function(x){x})
    s.obj@meta.data[, cluster.name] <- factor(s.obj@meta.data[, cluster.name],
                                                      levels=sort(help_sort_func(c_names)))
  }
  
  diff.pathways <- FindAllMarkers(s.obj, assay = "progeny", slot = "data", test.use = "wilcox") %>%
    subset(p_val_adj <= 0.05)
  row.names(diff.pathways) <- NULL
  
  s.obj@misc$progeny.diff.pathways <- diff.pathways
  
  if (save.RDS.s10 == TRUE){
    dir.create(file.path(path.to.output, "s10_output"), showWarnings = FALSE)
    saveRDS(object = cluster.markers, 
            file = file.path(path.to.output, "s10_output", 
                             paste0(PROJECT, ".PROGENy.s10.rds")))
  }
  return(s.obj)
}